#include "AlignmentReader.h"

namespace Mosaik {

	// constructor
	CAlignmentReader::CAlignmentReader(void)
		: mIsOpen(false)
		, mInStream(NULL)
		, mNumReads(0)
		, mCurrentRead(0)
		, mBuffer(NULL)
		, mBufferPtr(NULL)
		, mBufferLen(0)
		, mCompressionBuffer(NULL)
		, mCompressionBufferLen(0)
		, mPartitionSize(0)
		, mPartitionMembers(0)
		, mRefSeqLUT(NULL)
		, mStatus(AS_UNKNOWN)
		, mSeqTech(ST_UNKNOWN)
	{}

	// destructor
	CAlignmentReader::~CAlignmentReader(void) {
		if(mIsOpen)            Close();
		if(mBuffer)            delete [] mBuffer;
		if(mCompressionBuffer) delete [] mCompressionBuffer;

		// delete the reference sequence LUT
		for(unsigned short i = 0; i < mNumRefSeqs; i++) delete [] mRefSeqLUT[i];
		delete [] mRefSeqLUT;
	}

	// checks to see if this is truly an MOSAIK alignment archive
	bool CAlignmentReader::CheckFile(const string& filename, const bool showError) {

		// read in the first 6 characters
		char signature[7];
		signature[6] = 0;
		bool foundError = false;

		const char* MOSAIK_SIGNATURE = "MSKAA\2";

		// open the MOSAIK alignment archive
		FILE* checkStream = NULL;
		if(fopen_s(&checkStream, filename.c_str(), "rb") != 0) {
			if(showError) {
				printf("ERROR: Could not open %s when validating the alignment archive.\n", filename.c_str());
				exit(1);
			}

			foundError = true;
		}

		// retrieve the MOSAIK alignment archive signature
		if(!foundError) {

			// check if we were able to read 6 bytes
			if(fread(signature, 1, 6, checkStream) < 6) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) is not in the MOSAIK read format.\n", filename.c_str());
					exit(1);
				}

				foundError = true;
			}

			// check if the read signatures match
			if(!foundError && (strncmp(signature, MOSAIK_SIGNATURE, 5) != 0)) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) is not in the MOSAIK alignment format.\n", filename.c_str());
					exit(1);
				}

				foundError = true;
			}

			// check if the file format is from another version
			if(!foundError && (MOSAIK_SIGNATURE[5] != signature[5])) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) was created in another version of MosaikAligner. A new alignment archive is required.\n", filename.c_str());
					printf("       file version: %hu, expected version: %hu\n", signature[5], MOSAIK_SIGNATURE[5]);
					exit(1);
				}

				foundError = true;
			}
		}

		// close the file
		if(checkStream) fclose(checkStream);

		// return the appropriate values
		if(foundError) return false;
		return true;
	}

	// checks if the buffer is large enough to accomodate the requested size
	void CAlignmentReader::CheckBufferSize(unsigned char* &pBuffer, unsigned int& bufferLen, const unsigned int requestedBytes) {
		try {
			if(requestedBytes > bufferLen) {
				bufferLen = requestedBytes + 10;
				delete [] pBuffer;
				pBuffer = new unsigned char[bufferLen];
			}
		} catch(bad_alloc) {
			cout << "ERROR: Out of memory when allocating " << requestedBytes << endl;
			exit(1);
		}
	}

	// closes the alignment archive
	void CAlignmentReader::Close(void) {
		mIsOpen = false;
		fclose(mInStream);
	}

	// returns the number of bases in the archive
	uint64_t CAlignmentReader::GetNumBases(void) const {
		if(!mIsOpen) return 0;
		return mNumBases;
	}

	// returns the number of reads in the archive
	uint64_t CAlignmentReader::GetNumReads(void) const {
		if(!mIsOpen) return 0;
		return mNumReads;
	}

	// retrieves the read groups vector
	void CAlignmentReader::GetReadGroups(vector<ReadGroup>& readGroups) const {
		readGroups.resize(mReadGroups.size());
		vector<ReadGroup>::const_iterator rgIter;
		vector<ReadGroup>::iterator crgIter = readGroups.begin();

		// force a deep copy
		for(rgIter = mReadGroups.begin(); rgIter != mReadGroups.end(); rgIter++) {
			*crgIter = *rgIter;
		}
	}

	// retrieves the reference sequence data
	vector<ReferenceSequence>* CAlignmentReader::GetReferenceSequences(void) {
		return &mReferenceSequences;
	}

	// returns the sequencing technology
	SequencingTechnologies CAlignmentReader::GetSequencingTechnology(void) const {
		if(!mIsOpen) return ST_UNKNOWN;
		return mSeqTech;
	}
	
	// returns the alignment status
	AlignmentStatus CAlignmentReader::GetStatus(void) const {
		if(!mIsOpen) return AS_UNKNOWN;
		return mStatus;
	}

	// jumps to the block containing the specified reference index and position
	void CAlignmentReader::Jump(const unsigned int referenceIndex, const unsigned int referencePosition) {

		// ===============
		// parse the index
		// ===============

		if(mIndexOffset == 0) {
			cout << "ERROR: Cannot jump to the desired compressed block because the index offset was not set." << endl;
			exit(1);
		}

		fseek64(mInStream, mIndexOffset, SEEK_SET);

		// read the number of entries
		unsigned int numIndexEntries = 0;
		fread((char*)&numIndexEntries, SIZEOF_INT, 1, mInStream);

		// load the index
		const unsigned int requestedBytes = numIndexEntries * (SIZEOF_UINT64 + SIZEOF_INT + SIZEOF_SHORT);
		CheckBufferSize(mBuffer, mBufferLen, requestedBytes);
		fread(mBuffer, requestedBytes, 1, mInStream);

		// find the block containing the specified reference index and position
		unsigned int bufferOffset = 0;

		unsigned short index      = 0;
		unsigned int position     = 0;
		off_type offset           = 0;

		bool foundBlock = false;

		for(unsigned int i = 0; i < numIndexEntries; i++) {

			// retrieve the reference index
			memcpy((char*)&index, mBuffer + bufferOffset, SIZEOF_SHORT);
			bufferOffset += SIZEOF_SHORT;

			// store the reference position
			memcpy((char*)&position, mBuffer + bufferOffset, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;

			// store the file offset
			memcpy((char*)&offset, mBuffer + bufferOffset, SIZEOF_UINT64);
			bufferOffset += SIZEOF_UINT64;

			// keep going until we find a compression block that is past our desired index and position
			if(index > referenceIndex) foundBlock = true;
			if((index == referenceIndex) && (position >= referencePosition)) foundBlock = true;
			if(foundBlock) break;
		}

		if(!foundBlock) {
			cout << "ERROR: A suitable compression block was not found in the index." << endl;
			exit(1);
		}

		fseek64(mInStream, offset, SEEK_SET);
		mCurrentRead      = 0;
		mPartitionMembers = 0;
		mPartitionSize    = 0;
	}

	// loads the next read from the alignment archive
	bool CAlignmentReader::LoadNextRead(Mosaik::AlignedRead& ar) {

		if(!mIsOpen) {
			cout << "ERROR: An attempt was made to get reads from an alignment archive that hasn't been opened yet." << endl;
			exit(1);
		}

		// check if we have already processed all of the reads
		if(mCurrentRead >= mNumReads) return false;

		// ==================
		// read the partition
		// ==================

		if(mPartitionMembers == mPartitionSize) {

			// read the uncompressed partition entry size
			unsigned int uncompressedSize = 0;
			fread((char*)&uncompressedSize, SIZEOF_INT, 1, mInStream);

			if(feof(mInStream)) return false;

			// read the compressed partition entry size
			int compressedSize = 0;
			fread((char*)&compressedSize, SIZEOF_INT, 1, mInStream);

			// read the partition member size
			mPartitionMembers = 0;
			fread((char*)&mPartitionSize, SIZEOF_SHORT, 1, mInStream);

			// check the compression buffer size
			CheckBufferSize(mCompressionBuffer, mCompressionBufferLen, compressedSize);
			CheckBufferSize(mBuffer, mBufferLen, uncompressedSize);

			// read and uncompress the partition
			int numBytesRead = fread(mCompressionBuffer, 1, compressedSize, mInStream);

			if(numBytesRead != compressedSize) {
				cout << "ERROR: Tried to read " << compressedSize << " bytes, but received only " << numBytesRead << " bytes (" << mInputFilename << ") [read the partition: LoadNextRead]" << endl;
				exit(1);
			}

			int result = fastlz_decompress(mCompressionBuffer, compressedSize, mBuffer, mBufferLen);

			if(result == 0) {
				cout << "ERROR: Unable to properly uncompress the current data partition." << endl;
				exit(1);
			}

			// set the buffer pointer
			mBufferPtr = mBuffer;
		}

		// get the read name
		unsigned char readNameLen = *mBufferPtr;
		mBufferPtr++;

		ar.Name.resize(readNameLen);
		memcpy((void*)ar.Name.data(), mBufferPtr, readNameLen);
		mBufferPtr += readNameLen;

		// get the read group code
		memcpy((char*)&ar.ReadGroupCode, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;

		// get the status flag: long read, paired-end
		const unsigned status = *mBufferPtr;
		mBufferPtr++;

		const bool isLongRead  = ((status & 1) != 0 ? true : false);
		const bool isPairedEnd = ((status & 2) != 0 ? true : false);
		const bool hasMateInfo = ((status & 4) != 0 ? true : false);
		ar.IsLongRead = isLongRead;

		// get the number of mate 1 alignments
		unsigned int numMate1Alignments = 0;
		memcpy((char*)&numMate1Alignments, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;

		// get the number of mate 2 alignments
		unsigned int numMate2Alignments = 0;
		if(isPairedEnd) {
			memcpy((char*)&numMate2Alignments, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;
		}

		// =================================
		// deserialize each mate 1 alignment
		// =================================

		ar.Mate1Alignments.resize(numMate1Alignments);
		ReadAlignments(ar.Mate1Alignments, isLongRead, hasMateInfo);

		// =================================
		// deserialize each mate 2 alignment
		// =================================

		ar.Mate2Alignments.resize(numMate2Alignments);
		ReadAlignments(ar.Mate2Alignments, isLongRead, hasMateInfo);

		// increment the read counter
		mCurrentRead++;
		mPartitionMembers++;

		return true;
	}

	// deserializes each alignment and stores them in the supplied vector
	void CAlignmentReader::ReadAlignments(vector<Alignment>& alignments, const bool isLongRead, const bool hasMateInfo) {
		vector<Alignment>::iterator alIter;
		for(alIter = alignments.begin(); alIter != alignments.end(); alIter++) ReadAlignment(*alIter, isLongRead, hasMateInfo);
	}

	// deserialize the alignment
	void CAlignmentReader::ReadAlignment(Alignment& al, const bool isLongRead, const bool hasMateInfo) {

		// get the reference sequence start position
		memcpy((char*)&al.ReferenceBegin, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;

		// get the reference sequence end position
		memcpy((char*)&al.ReferenceEnd, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;

		// get the reference sequence index
		memcpy((char*)&al.ReferenceIndex, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;
		al.ReferenceName = mRefSeqLUT[al.ReferenceIndex];

		// get the alignment quality
		al.Quality = *mBufferPtr;
		mBufferPtr++;

		// skip the alternate alignment quality
		al.AlternateQuality = *mBufferPtr;
		mBufferPtr++;

		// get the read orientation
		const unsigned char readOrientation = *mBufferPtr;
		mBufferPtr++;

		al.IsReverseComplement     = false;
		al.IsMateReverseComplement = false;
		if((readOrientation & 1) == 1) al.IsReverseComplement     = true;
		if((readOrientation & 2) == 1) al.IsMateReverseComplement = true;

		// get the mate reference info
		if(hasMateInfo) {

			// get the mate reference sequence start position
			memcpy((char*)&al.MateReferenceBegin, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;

			// get the mate reference sequence end position
			memcpy((char*)&al.MateReferenceEnd, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;

			// get the mate reference sequence index
			memcpy((char*)&al.MateReferenceIndex, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;
		}

		unsigned short pairwiseLength = 0;

		if(isLongRead) {

			// get the pairwise length
			memcpy((char*)&pairwiseLength, mBufferPtr, SIZEOF_SHORT);
			mBufferPtr += SIZEOF_SHORT;

			// get the query begin
			memcpy((char*)&al.QueryBegin, mBufferPtr, SIZEOF_SHORT);
			mBufferPtr += SIZEOF_SHORT;

			// get the query end
			memcpy((char*)&al.QueryEnd, mBufferPtr, SIZEOF_SHORT);
			mBufferPtr += SIZEOF_SHORT;

		} else {

			// get the pairwise length
			pairwiseLength = *mBufferPtr;
			mBufferPtr++;

			// get the query begin
			al.QueryBegin = *mBufferPtr;
			mBufferPtr++;

			// get the query end
			al.QueryEnd = *mBufferPtr;
			mBufferPtr++;
		}

		// get the pairwise reference bases
		al.Reference.resize(pairwiseLength);
		memcpy((void*)al.Reference.data(), mBufferPtr, pairwiseLength);
		mBufferPtr += pairwiseLength;

		// get the pairwise query bases
		al.Query.resize(pairwiseLength);
		memcpy((void*)al.Query.data(), mBufferPtr, pairwiseLength);
		mBufferPtr += pairwiseLength;

		// get the pairwise query base qualities
		const unsigned short bqLength = al.QueryEnd - al.QueryBegin + 1;
		al.BaseQualities.resize(bqLength);
		memcpy((void*)al.BaseQualities.data(), mBufferPtr, bqLength);
		mBufferPtr += bqLength;

		// DEBUG
		//cout << "reference index: " << al.ReferenceIndex << ", begin: " << al.ReferenceBegin << ", end: " << al.ReferenceEnd << endl;
		//cout << "query begin: " << al.QueryBegin << ", end: " << al.QueryEnd << ", pairwise length: " << pairwiseLength << endl << endl;
	}

	// opens the alignment archive
	void CAlignmentReader::Open(const string& filename) {

		if(mIsOpen) {
			cout << "ERROR: An attempt was made to open an already open alignment archive." << endl;
			exit(1);
		}

		mInputFilename = filename;

		if(fopen_s(&mInStream, filename.c_str(), "rb") != 0) {
			cout << "ERROR: Could not open the compressed alignment archive (" << mInputFilename << ") for reading." << endl;
			exit(1);
		}

		mIsOpen = true;

		// ===============
		// read the header
		// ===============

		// MOSAIK_SIGNATURE[6]	   0  -  5
		// STATUS[1]               6  -  6
		// SEQUENCE_TECHNOLOGY[2]  7  -  8
		// ARCHIVE_DATE[8]		   9  - 16
		// NUM_REFERENCE_SEQS[4]   17 - 20
		// NUM_READ_GROUPS[4]      21 - 24
		// NUM_READS[8]            25 - 32
		// NUM_BASES[8]            33 - 40
		// REFERENCES_OFFSET[8]    41 - 48
		// REFERENCE_GAP_OFFSET[8] 49 - 57
		// INDEX_OFFSET[8]         58 - 63
		// RESERVED[8]             64 - 71
		// READ_GROUPS[*]

		// skip the MOSAIK signature
		const unsigned char SIGNATURE_LENGTH = 6;
		fseek64(mInStream, SIGNATURE_LENGTH, SEEK_SET);

		// retrieve the alignment file status
		mStatus = (AlignmentStatus)fgetc(mInStream);

		// retrieve the sequencing technology
		fread((char*)&mSeqTech, SIZEOF_SHORT, 1, mInStream);

		// skip the archive date
		fseek64(mInStream, SIZEOF_UINT64, SEEK_CUR);

		// retrieve the number of reference sequences
		fread((char*)&mNumRefSeqs, SIZEOF_INT, 1, mInStream);

		// retrieve the number of read groups
		unsigned int numReadGroups;
		fread((char*)&numReadGroups, SIZEOF_INT, 1, mInStream);

		// retrieve the number of reads
		fread((char*)&mNumReads, SIZEOF_UINT64, 1, mInStream);

		// retrieve the number of bases
		fread((char*)&mNumBases, SIZEOF_UINT64, 1, mInStream);

		// retrieve the references offset
		off_type referencesOffset = 0;
		fread((char*)&referencesOffset, SIZEOF_UINT64, 1, mInStream);

		// retrieve the reference gaps offset
		fread((char*)&mReferenceGapOffset, SIZEOF_UINT64, 1, mInStream);

		// retrieve the index offset
		fread((char*)&mIndexOffset, SIZEOF_UINT64, 1, mInStream);

		// skip the reserved
		fseek64(mInStream, SIZEOF_UINT64, SEEK_CUR);

		//// DEBUG
		//cout << "mStatus:             " << (short)mStatus << endl;
		//cout << "mSeqTech:            " << mSeqTech << endl;
		//cout << "mNumRefSeqs:         " << mNumRefSeqs << endl;
		//cout << "numReadGroups:       " << numReadGroups << endl;
		//cout << "mNumReads:           " << mNumReads << endl;
		//cout << "mNumBases:           " << mNumBases << endl;
		//cout << "referencesOffset:    " << referencesOffset << endl;
		//cout << "mReferenceGapOffset: " << mReferenceGapOffset << endl;
		//cout << "mIndexOffset:        " << mIndexOffset << endl;
		//cout << "reserved:            " << reserved << endl;
		//exit(1);

		// retrieve the read groups
		mReadGroups.resize(numReadGroups);

		vector<ReadGroup>::iterator rgIter;
		for(rgIter = mReadGroups.begin(); rgIter != mReadGroups.end(); rgIter++) {

			// read the metadata string lengths
			const unsigned char centerNameLen   = (unsigned char)fgetc(mInStream);
			const unsigned char libraryNameLen  = (unsigned char)fgetc(mInStream);
			const unsigned char platformUnitLen = (unsigned char)fgetc(mInStream);
			const unsigned char readGroupIDLen  = (unsigned char)fgetc(mInStream);
			const unsigned char sampleNameLen   = (unsigned char)fgetc(mInStream);

			unsigned short descriptionLen = 0;
			fread((char*)&descriptionLen, SIZEOF_SHORT, 1, mInStream);
			fread((char*)&rgIter->SequencingTechnology, SIZEOF_SHORT, 1, mInStream);
			fread((char*)&rgIter->MedianFragmentLength, SIZEOF_INT, 1, mInStream);

			// skip the reserved bytes
			fseek64(mInStream, SIZEOF_UINT64, SEEK_CUR);

			rgIter->CenterName.resize(centerNameLen);
			rgIter->LibraryName.resize(libraryNameLen);
			rgIter->PlatformUnit.resize(platformUnitLen);
			rgIter->ReadGroupID.resize(readGroupIDLen);
			rgIter->SampleName.resize(sampleNameLen);
			rgIter->Description.resize(descriptionLen);

			// read the metadata strings
			fread((void*)rgIter->CenterName.data(),   centerNameLen,   1, mInStream);
			fread((void*)rgIter->Description.data(),  descriptionLen,  1, mInStream);
			fread((void*)rgIter->LibraryName.data(),  libraryNameLen,  1, mInStream);
			fread((void*)rgIter->PlatformUnit.data(), platformUnitLen, 1, mInStream);
			fread((void*)rgIter->ReadGroupID.data(),  readGroupIDLen,  1, mInStream);
			fread((void*)rgIter->SampleName.data(),   sampleNameLen,   1, mInStream);

			//// DEBUG
			//cout << "center name:            " << rgIter->CenterName << endl;
			//cout << "description:            " << rgIter->Description << endl;
			//cout << "library name:           " << rgIter->LibraryName << endl;
			//cout << "platform unit:          " << rgIter->PlatformUnit << endl;
			//cout << "read group ID:          " << rgIter->ReadGroupID << endl;
			//cout << "sample name:            " << rgIter->SampleName << endl;
			//cout << "sequencing technology:  " << rgIter->SequencingTechnology << endl;
			//cout << "median fragment length: " << rgIter->MedianFragmentLength << endl;
		}

		// store the reads offset
		mReadsOffset = ftell64(mInStream);

		// ============================
		// read the reference sequences
		// ============================

		// jump to the reference sequence section
		fseek64(mInStream, referencesOffset, SEEK_SET);

		mReferenceSequences.resize(mNumRefSeqs);
		mRefSeqLUT = new char*[mNumRefSeqs];

		unsigned int currentRefSeq = 0;
		vector<ReferenceSequence>::iterator rsIter;
		for(rsIter = mReferenceSequences.begin(); rsIter != mReferenceSequences.end(); rsIter++, currentRefSeq++) {

			// REFERENCE_SEQ_NAME_LEN[1]                0 -  0 
			// REFERENCE_SEQ_SPECIES_LEN[1]             1 -  1
			// REFERENCE_SEQ_GENOME_ASSEMBLY_ID_LEN[1]  2 -  2
			// REFERENCE_SEQ_URI_LEN[1]                 3 -  3
			// REFERENCE_SEQ_NUM_BASES[4]               4 -  7
			// REFERENCE_SEQ_SEQ_OFFSET[8]              8 - 15
			// REFERENCE_SEQ_MD5[16]                   16 - 31
			// REFERENCE_SEQ_NAME[X]                   32 - XX
			// REFERENCE_SEQ_SPECIES[X]
			// REFERENCE_SEQ_GENOME_ASSEMBLY_ID[X]
			// REFERENCE_SEQ_URI[X]

			// read the name length
			const unsigned char nameLen = fgetc(mInStream);

			// read the species length
			const unsigned char speciesLen = fgetc(mInStream);

			// read the genome assembly id length
			const unsigned char genomeAssemblyIDLen = fgetc(mInStream);

			// read the uri length
			const unsigned char uriLen = fgetc(mInStream);

			// read the number of bases
			fread((char*)&rsIter->NumBases, SIZEOF_INT, 1, mInStream);

			// write the number of aligned reads
			fread((char*)&rsIter->NumAligned, SIZEOF_UINT64, 1, mInStream);

			// read the MD5 checksum
			rsIter->MD5.resize(32);
			char* pBuffer = (char*)rsIter->MD5.data();
			fread(pBuffer, 32, 1, mInStream);

			// read the reference name
			rsIter->Name.resize(nameLen);
			pBuffer = (char*)rsIter->Name.data();
			fread(pBuffer, nameLen, 1, mInStream);

			mRefSeqLUT[currentRefSeq] = new char[nameLen + 1];
			memcpy(mRefSeqLUT[currentRefSeq], pBuffer, nameLen);
			mRefSeqLUT[currentRefSeq][nameLen] = 0;

			// read the species name
			if(speciesLen > 0) {
				rsIter->Species.resize(speciesLen);
				pBuffer = (char*)rsIter->Species.data();
				fread(pBuffer, speciesLen, 1, mInStream);
			}

			// read the genome assembly ID
			if(genomeAssemblyIDLen > 0) {
				rsIter->GenomeAssemblyID.resize(genomeAssemblyIDLen);
				pBuffer = (char*)rsIter->GenomeAssemblyID.data();
				fread(pBuffer, genomeAssemblyIDLen, 1, mInStream);
			}

			// read the URI
			if(uriLen > 0) {
				rsIter->URI.resize(uriLen);
				pBuffer = (char*)rsIter->URI.data();
				fread(pBuffer, uriLen, 1, mInStream);
			}

			//// DEBUG
			//cout << "# bases:                " << rsIter->NumBases << endl;
			//cout << "md5:                    " << rsIter->MD5 << endl;
			//cout << "name:                   " << rsIter->Name << endl;
			//cout << "species:                " << rsIter->Species << endl;
			//cout << "genome assembly ID:     " << rsIter->GenomeAssemblyID << endl;
			//cout << "URI:                    " << rsIter->URI << endl;
		}

		// restore our file position
		Rewind();
	}

	// sets the file pointer to the beginning of the read data
	void CAlignmentReader::Rewind(void) {
		fseek64(mInStream, mReadsOffset, SEEK_SET);
		mCurrentRead      = 0;
		mPartitionMembers = 0;
		mPartitionSize    = 0;
	}
}
