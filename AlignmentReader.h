#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ALIGNMENT_NO_MATE_INFO 0xffffffff

// define some data type sizes
#define SIZEOF_CHAR          1
#define SIZEOF_SHORT         2
#define SIZEOF_INT           4
#define SIZEOF_UINT64        8

// define our 64-bit data type
#ifdef WIN32
typedef unsigned long long uint64_t;
#endif

typedef unsigned short SequencingTechnologies;
typedef unsigned char AlignmentStatus;

#define ST_UNKNOWN               0
#define ST_454                   1
#define ST_HELICOS               2
#define ST_ILLUMINA              4
#define ST_PACIFIC_BIOSCIENCES   8
#define ST_SOLID                16
#define ST_SANGER               32

#define AS_UNKNOWN           0 
#define AS_SINGLE_END_READ   1 // transferred from the read format
#define AS_PAIRED_END_READ   2 // transferred from the read format
#define AS_UNSORTED_READ     4 // expected in MosaikAligner data
#define AS_SORTED_ALIGNMENT  8 // expected in MosaikSort data
#define AS_ALL_MODE         16 // enables non-unique PE resolution
#define AS_UNIQUE_MODE      32 // disables non-unique PE resolution

// add some large file support
#ifdef WIN32
#define ftell64(a)     _ftelli64(a)
#define fseek64(a,b,c) _fseeki64(a,b,c)
typedef __int64 off_type;
#elif SPARC
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef off64_t off_type;
#elif MACOSX
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef fpos_t off_type;
#else
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef off_t off_type;
#endif

// add my implementation of the "safe" functions
#ifndef WIN32
typedef int errno_t;
#define EINVAL 22
inline errno_t fopen_s(FILE** pFile, const char *filename, const char *mode) {
	*pFile = fopen(filename, mode);
	if(*pFile) return 0;
	return EINVAL;
}
#endif

using namespace std;

namespace Mosaik {

	struct Alignment {
		unsigned int MateReferenceBegin;
		unsigned int MateReferenceEnd;
		unsigned int MateReferenceIndex;
		unsigned int ReferenceBegin;
		unsigned int ReferenceEnd;
		unsigned short ReferenceIndex;
		unsigned short QueryBegin;
		unsigned short QueryEnd;
		unsigned char AlternateQuality;
		unsigned char Quality;
		bool IsReverseComplement;
		bool IsMateReverseComplement;
		char* ReferenceName;
		string Reference;
		string Query;
		string BaseQualities;

		// constructor
		Alignment(void)
			: MateReferenceIndex(ALIGNMENT_NO_MATE_INFO)
			, ReferenceIndex(0)
			, Quality(0)
			, IsReverseComplement(false)
			, IsMateReverseComplement(false)
			, ReferenceName(NULL)
		{}
	};

	struct AlignedRead {
		unsigned int ReadGroupCode;
		string Name;
		vector<Alignment> Mate1Alignments;
		vector<Alignment> Mate2Alignments;
		bool IsLongRead;

		// constructor
		AlignedRead()
			: ReadGroupCode(0)
			, IsLongRead(false)
		{}
	};

	struct ReferenceSequence {
		off_type BasesOffset;
		uint64_t NumAligned;
		unsigned int Begin;
		unsigned int End;
		unsigned int NumBases;
		string Name;
		string Bases;
		string GenomeAssemblyID;
		string Species;
		string MD5;
		string URI;

		// constructor
		ReferenceSequence()
			: BasesOffset(0)
			, NumAligned(0)
			, Begin(0)
			, End(0)
			, NumBases(0)
		{}
	};

	struct ReadGroup {
		unsigned int MedianFragmentLength;
		unsigned int ReadGroupCode;
		SequencingTechnologies SequencingTechnology;
		string CenterName;
		string Description;
		string LibraryName;
		string PlatformUnit;
		string ReadGroupID;
		string SampleName;

		// constructor
		ReadGroup(void)
			: MedianFragmentLength(0)
			, ReadGroupCode(0)
			, SequencingTechnology(ST_UNKNOWN)
		{}
	};

	class CAlignmentReader {
	public:
		// constructor
		CAlignmentReader(void);
		// destructor
		~CAlignmentReader(void);
		// checks if the buffer is large enough to accomodate the requested size
		static void CheckBufferSize(unsigned char* &pBuffer, unsigned int& bufferLen, const unsigned int requestedBytes);
		// checks to see if this is truly a MOSAIK alignment archive
		static bool CheckFile(const string& filename, const bool showError);
		// closes the alignment archive
		void Close(void);
		// returns the number of bases in the archive
		uint64_t GetNumBases(void) const;
		// returns the number of reads in the archive
		uint64_t GetNumReads(void) const;
		// retrieves the read groups vector
		void GetReadGroups(vector<ReadGroup>& readGroups) const;
		// retrieves the reference sequence data
		vector<ReferenceSequence>* GetReferenceSequences(void);
		// gets the alignment archive sequencing technology
		SequencingTechnologies GetSequencingTechnology(void) const;
		// retrieves the file status
		AlignmentStatus GetStatus(void) const;
		// jumps to the block containing the specified reference index and position
		void Jump(const unsigned int referenceIndex, const unsigned int referencePosition);
		// loads the next read from the alignment archive
		bool LoadNextRead(Mosaik::AlignedRead& ar);
		// opens the alignment archive
		void Open(const string& filename);
		// sets the file pointer to the beginning of the read data
		void Rewind(void);

	private:
		// deserializes each alignment and stores them in the supplied vector
		void ReadAlignments(vector<Alignment>& alignments, const bool isLongRead, const bool hasMateInfo);
		// deserialize the alignment
		void ReadAlignment(Alignment& al, const bool isLongRead, const bool hasMateInfo);
		// denotes the status of the output stream
		bool mIsOpen;
		// our compressed output stream
		FILE* mInStream;
		// stores the archive read count
		uint64_t mNumReads;
		uint64_t mNumBases;
		// stores the current read number
		uint64_t mCurrentRead;
		// stores the file offsets
		off_type mReadsOffset;
		off_type mReferenceGapOffset;
		off_type mIndexOffset;
		// our input buffer
		unsigned char* mBuffer;
		unsigned char* mBufferPtr;
		unsigned int mBufferLen;
		// our input compression buffer
		unsigned char* mCompressionBuffer;
		unsigned int mCompressionBufferLen;
		// our input filename
		string mInputFilename;
		// our partitioning setup
		unsigned short mPartitionSize;
		unsigned short mPartitionMembers;
		// our reference sequence LUT
		char** mRefSeqLUT;
		unsigned short mNumRefSeqs;
		// our reference sequences
		vector<ReferenceSequence> mReferenceSequences;
		// our read groups
		vector<ReadGroup> mReadGroups;
		// our file status
		AlignmentStatus mStatus;
		// our sequencing technology
		SequencingTechnologies mSeqTech;
	};
}

// add a prototype for the Ariya Hidayat's FastLZ library
extern "C" {
	int fastlz_decompress(const void* input, int length, void* output, int maxout);
}
