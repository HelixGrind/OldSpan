#include "AlignedPositionsReader.h"

namespace MosaikReadFormat {

	CAlignedPositionsReader::CAlignedPositionsReader(const char* filename)
		: mIsStreamOpen(false)
		, mInStream(NULL)
		, mBuffer(NULL)
		, mBufferLen(256)
		, mReadBuffer(NULL)
		, mReadBufferLen(256)
		, mCurrentReadNum(0)
	{
		// initialize the read and index buffer
		try {
			mBuffer = new unsigned char[mBufferLen];
		} catch(bad_alloc) {
			double requestedBytes = SIZEOF_CHAR * mBufferLen;
			cout << "ERROR: Out of memory when allocating " << requestedBytes << " for the Mosaik aligned read format buffer." << endl;
			exit(1);
		}

		try {
			mReadBuffer = new char[mReadBufferLen];
		} catch(bad_alloc) {
			double requestedBytes = SIZEOF_CHAR * mReadBufferLen;
			cout << "ERROR: Out of memory when allocating " << requestedBytes << " for the Mosaik aligned read format buffer." << endl;
			exit(1);
		}

		// open the file
		errno_t err;
		if((err = fopen_s(&mInStream, filename, "rb")) != 0) {
			cout << "ERROR: Could not open " << filename << " (Mosaik read format) for reading." << endl;
			exit(1);
		}
		mIsStreamOpen = true;

		// parse the header
		ParseHeader();
	}

	CAlignedPositionsReader::~CAlignedPositionsReader(void) {
		if(mIsStreamOpen) Close();
	}

	// checks to see if this is truly a Mosaik aligned read report
	bool CAlignedPositionsReader::CheckMosaikAlignedReadArchive(const char* filename, bool exitOnError) {
		
		// read in the first four characters
		char signature[7];
		signature[6] = 0;

		char* MOSAIK_SIGNATURE = "MSKAR\0";

		// open the Mosaik aligned read report
		FILE* checkStream = NULL;
		errno_t err;

		if((err = fopen_s(&checkStream, filename, "rb")) != 0) {
			if(exitOnError) {
				cout << "ERROR: Could not open " << filename << " when checking the aligned read report signature." << endl;
				exit(1);
			}
			return false;
		}

		// retrieve the Mosaik aligned read report signature
		fread(signature, 6, 1, checkStream);

		// check if the read signatures match
		if(strncmp(signature, MOSAIK_SIGNATURE, 6) != 0) {
			if(exitOnError) {
				cout << "ERROR: It seems that the input file (" << filename << ") is not in the Mosaik aligned read report format." << endl;
				exit(1);
			}
			return false;
		}

		// close the file
		fclose(checkStream);

		return true;
	}
  
  // returns the number of reads in this file
  mosaik_64u_t CAlignedPositionsReader::GetNumReads(void) const {
    return mStatistics.NumReads;
  }

	// closes the aligned positions file
	void CAlignedPositionsReader::Close(void) {
		if(mIsStreamOpen) fclose(mInStream);
		mIsStreamOpen      = false;
	}

	// parses the aligned positions file header
	void CAlignedPositionsReader::ParseHeader(void) {

		// initialize
		unsigned int bufferOffset = 6;

		// read the header from disk
		size_t numBytesRead = fread(mBuffer, 1, ALIGNED_READ_ARCHIVE_HEADER_BYTES, mInStream);

		if(numBytesRead != ALIGNED_READ_ARCHIVE_HEADER_BYTES) {
			cout << "ERROR: Tried to read " << ALIGNED_READ_ARCHIVE_HEADER_BYTES << " bytes from the aligned positions file, but only read " << numBytesRead << " bytes." << endl;
			exit(1);
		}

		// extract the anchors offset
		memcpy((char*)&mStatistics.AnchorsOffset, mBuffer + bufferOffset, SIZEOF_MOSAIK_64);
		bufferOffset += SIZEOF_MOSAIK_64;

		//// extract the read names offset
		//memcpy((char*)&mStatistics.ReadNamesOffset, mBuffer + bufferOffset, SIZEOF_MOSAIK_64);
		//bufferOffset += SIZEOF_MOSAIK_64;
		//mStatistics.CurrentReadNameOffset = mStatistics.ReadNamesOffset;

		// extract the number of anchors
		memcpy((char*)&mStatistics.NumAnchors, mBuffer + bufferOffset, SIZEOF_INT);
		bufferOffset += SIZEOF_INT;

		// extract the number of reads
		memcpy((char*)&mStatistics.NumReads, mBuffer + bufferOffset, SIZEOF_MOSAIK_64);
		bufferOffset += SIZEOF_MOSAIK_64;

		// extract the archive timestamp
		memcpy((char*)&mStatistics.Timestamp, mBuffer + bufferOffset, SIZEOF_MOSAIK_64);

		//char* MOSAIK_SIGNATURE = "MSKAP\0";

		// DEBUG
		//cout << "# Reads:           " << mStatistics.NumReads << endl;
		//cout << "# Anchors:         " << mStatistics.NumAnchors << endl;
		//cout << "Anchors offset:    " << mStatistics.AnchorsOffset << endl;
		//cout << "Timestamp:         " << mStatistics.Timestamp << endl;
		//exit(1);
	}

	// checks if the buffer size is large enough to accomodate the requested size
	void CAlignedPositionsReader::CheckBufferSize(const unsigned int requestedBytes) {
		try {
			if(requestedBytes > mBufferLen) {
				mBufferLen = requestedBytes + 10;
				delete [] mBuffer;
				mBuffer = new unsigned char[mBufferLen];
			}
		} catch(bad_alloc) {
			cout << "ERROR: Out of memory when allocating " << requestedBytes << " for the Mosaik aligned read format buffer." << endl;
			exit(1);
		}
	}

	// checks if the read buffer size is large enough to accomodate the requested size
	void CAlignedPositionsReader::CheckReadBufferSize(const unsigned int requestedBytes) {
		try {
			if(requestedBytes > mReadBufferLen) {
				mReadBufferLen = requestedBytes + 10;
				delete [] mReadBuffer;
				mReadBuffer = new char[mReadBufferLen];
			}
		} catch(bad_alloc) {
			cout << "ERROR: Out of memory when allocating " << requestedBytes << " for the Mosaik aligned read format buffer." << endl;
			exit(1);
		}
	}

	// loads the next read from the aligned positions file
	bool CAlignedPositionsReader::LoadNextRead(vector<AlignedPosition>& positions, unsigned short& readLength, string& readName) {

		// have we already read our reads?
		if(mCurrentReadNum >= mStatistics.NumReads) return false;

		// clear our vector
		positions.clear();

		// get the number of aligned fragments & read length
		size_t numBytesRead = fread(mBuffer, 1, 10, mInStream);

		if(numBytesRead != 10) {
			cout << "ERROR: Tried to read " << 10 << " bytes from the aligned positions file, but only read " << numBytesRead << " bytes." << endl;
			exit(1);
		}

		// extract the entry size
		unsigned int entrySize = 0;
		memcpy((char*)&entrySize, mBuffer, SIZEOF_INT);
		unsigned int bufferOffset = SIZEOF_INT;

		// extract the number of aligned positions
		unsigned int numAlignedPositions = 0;
		memcpy((char*)&numAlignedPositions, mBuffer + bufferOffset, SIZEOF_INT);
		bufferOffset += SIZEOF_INT;
		
		// extract the read length
		memcpy((char*)&readLength, mBuffer + bufferOffset, SIZEOF_SHORT);
		bufferOffset += SIZEOF_SHORT;

		// DEBUG
		//cout << "entry size:          " << entrySize << endl;
		//cout << "# aligned positions: " << numAlignedPositions << endl;
		//cout << "read length:         " << readLength << endl;
		//exit(1);

		// make sure our buffer size is large enough
		unsigned int entryBytesLeft = entrySize - 10;
		CheckBufferSize(entryBytesLeft);

		// get the rest of the entry data
		numBytesRead = fread(mBuffer, 1, entryBytesLeft, mInStream);
		
		if(numBytesRead != entryBytesLeft) {
			cout << "ERROR: Tried to read " << entryBytesLeft << " bytes from the aligned positions file, but only read " << numBytesRead << " bytes." << endl;
			exit(1);
		}
	
		// extract the read name length
		bufferOffset = 0;
		unsigned char readNameLength = mBuffer[bufferOffset++];

		// extract the read name
		CheckReadBufferSize(readNameLength);
		memcpy(mReadBuffer, mBuffer + bufferOffset, readNameLength);
		bufferOffset += readNameLength;
		mReadBuffer[readNameLength] = 0;
		readName = mReadBuffer;

		// ====================================
		// extract all of the aligned positions 
		// ====================================

		for(unsigned int i = 0; i < numAlignedPositions; i++) {

			AlignedPosition arp;

			// extract the anchor start
			memcpy((char*)&arp.AnchorBegin, mBuffer + bufferOffset, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;

			// extract the anchor end
			memcpy((char*)&arp.AnchorEnd, mBuffer + bufferOffset, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;

			// extract the query start
			memcpy((char*)&arp.QueryBegin, mBuffer + bufferOffset, SIZEOF_SHORT);
			bufferOffset += SIZEOF_SHORT;

			// extract the query end
			memcpy((char*)&arp.QueryEnd, mBuffer + bufferOffset, SIZEOF_SHORT);
			bufferOffset += SIZEOF_SHORT;

			// set the reverse complement flag
			if(mBuffer[bufferOffset] == 1) arp.IsReverseComplement = true;
				else arp.IsReverseComplement = false;
			bufferOffset++;

			// extract the number of substitutions
			memcpy((char*)&arp.NumSubstitutions, mBuffer + bufferOffset, SIZEOF_SHORT);
			bufferOffset += SIZEOF_SHORT;

			// extract the number of insertions
			memcpy((char*)&arp.NumInsertions, mBuffer + bufferOffset, SIZEOF_SHORT);
			bufferOffset += SIZEOF_SHORT;

			// extract the number of deletions
			memcpy((char*)&arp.NumDeletions, mBuffer + bufferOffset, SIZEOF_SHORT);
			bufferOffset += SIZEOF_SHORT;

			// extract the pairwise length
			unsigned short pairwiseLength = 0;
			memcpy((char*)&pairwiseLength, mBuffer + bufferOffset, SIZEOF_SHORT);
			bufferOffset += SIZEOF_SHORT;

			// extract the anchor
			CheckReadBufferSize(pairwiseLength);
			memcpy(mReadBuffer, mBuffer + bufferOffset, pairwiseLength);
			bufferOffset += pairwiseLength;
			mReadBuffer[pairwiseLength] = 0;
			arp.Anchor = mReadBuffer;

			// extract the query
			memcpy(mReadBuffer, mBuffer + bufferOffset, pairwiseLength);
			bufferOffset += pairwiseLength;
			mReadBuffer[pairwiseLength] = 0;
			arp.Query = mReadBuffer;

			// add the structure to the vector
			positions.push_back(arp);

			// DEBUG
			//cout << "anchor:             " << arp.Anchor << endl;
			//cout << "query:              " << arp.Query << endl;
			//cout << "anchor begin:       " << arp.AnchorBegin << endl;
			//cout << "anchor end:         " << arp.AnchorEnd << endl;
			//cout << "query begin:        " << arp.QueryBegin << endl;
			//cout << "query end:          " << arp.QueryEnd << endl;
			//cout << "# of substitutions: " << arp.NumSubstitutions << endl;
			//cout << "# of insertions:    " << arp.NumInsertions << endl;
			//cout << "# of deletions:     " << arp.NumDeletions << endl;
			//cout << "orientation:        " << (arp.IsReverseComplement ? "reverse" : "forward") << endl;
			//exit(1);
		}
		//exit(1);

		// increment our read counter
		mCurrentReadNum++;

		//// store the current offset and jumpt to the read name offset
		//off_type previousOffset = ftell64(mInStream);
		//fseek64(mInStream, mStatistics.CurrentReadNameOffset, SEEK_SET);

		//// retrieve the read name
		//unsigned char readNameLength = 0;
		//fread((char*)&readNameLength, 1, 1, mInStream);
		//fread(mBuffer, 1, readNameLength, mInStream);
		//mBuffer[readNameLength] = 0;
		//readName = (char*)mBuffer;

		//// store the current read name offset and restore previous offset
		//mStatistics.CurrentReadNameOffset = ftell64(mInStream);
		//fseek64(mInStream, previousOffset, SEEK_SET);

		// returns
		return true;
	}

	// loads the anchors for the aligned positions file
	void CAlignedPositionsReader::LoadAnchors(vector<Anchor>& anchors) {

		// grab the current offset so we can restore it later
		off_type previousOffset = ftell64(mInStream);

		// jump to the anchors offset
		fseek64(mInStream, mStatistics.AnchorsOffset, SEEK_SET);

		// clear our vector
		anchors.clear();

		// =======================
		// copy all of the anchors
		// =======================

		for(unsigned int i = 0; i < mStatistics.NumAnchors; i++) {

			// define our variables
			Anchor a;
			unsigned char anchorNameLength = 0;
			unsigned int bufferOffset      = 0;

			// get the anchor name length
			fread((char*)&anchorNameLength, 1, 1, mInStream);

			// calculate how many bytes we need
			unsigned int numBytes = (ALIGNED_READ_ARCHIVE_ANCHOR_BYTES - 1) + anchorNameLength;

			// retrieve the anchor from disk
			size_t numBytesRead = fread(mBuffer, 1, numBytes, mInStream);
		
			if(numBytesRead != numBytes) {
				cout << "ERROR: Tried to read " << numBytes << " bytes from the aligned positions file, but only read " << numBytesRead << " bytes." << endl;
				exit(1);
			}

			// DEBUG
			//for(unsigned int j = bufferOffset; j < numBytes; j++) {
			//	cout << std::hex << (int)mBuffer[j] << " ";
			//}
			//exit(1);

			// extract the anchor length
			bufferOffset = anchorNameLength;
			memcpy((char*)&a.Length, mBuffer + bufferOffset, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;

			// extract the anchor start position
			memcpy((char*)&a.Begin, mBuffer + bufferOffset, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;

			// extract the anchor end position
			memcpy((char*)&a.End, mBuffer + bufferOffset, SIZEOF_INT);
			
			// extract the anchor name
			mBuffer[anchorNameLength] = 0;
			a.Name = (char*)mBuffer;

			// DEBUG
			//cout << "name:   " << a.Name << endl;
			//cout << "length: " << std::hex << a.Length << endl;
			//cout << "begin:  " << a.Begin << endl;
			//cout << "end:    " << a.End << endl;
			//exit(1);

			// add the structure to the vector
			anchors.push_back(a);
		}

		// restore the previous offset
		fseek64(mInStream, previousOffset, SEEK_SET);
	}
}
