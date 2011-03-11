#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <time.h>
#include "LargeFileSupport.h"
#include "SafeFunctions.h"

using namespace std;

// NOTE: it's a bit redundant, but let's make this class as standalone
// as possible to allow people to easily use this in other projects.
typedef unsigned long long mosaik_64u_t;

#define SIZEOF_CHAR      1
#define SIZEOF_SHORT     2
#define SIZEOF_INT       4
#define SIZEOF_MOSAIK_64 8

#define ALIGNED_READ_ARCHIVE_HEADER_BYTES 34
#define ALIGNED_READ_ARCHIVE_READ_BYTES   21
#define ALIGNED_READ_ARCHIVE_ANCHOR_BYTES 13

namespace MosaikReadFormat {

	class CAlignedPositionsReader {
	public:
		CAlignedPositionsReader(const char* filename);
		~CAlignedPositionsReader(void);
		// define our aligned read position structure
		struct AlignedPosition {
			unsigned int AnchorBegin;
			unsigned int AnchorEnd;
			unsigned short QueryBegin;
			unsigned short QueryEnd;
			unsigned short NumSubstitutions;
			unsigned short NumDeletions;
			unsigned short NumInsertions;
			bool IsReverseComplement;
			string AnchorName;
			string Anchor;
			string Query;

			bool operator<(const AlignedPosition& ar) const {
				if(AnchorBegin != ar.AnchorBegin) return AnchorBegin < ar.AnchorBegin;
				if(QueryBegin  != ar.QueryBegin)  return QueryBegin  < ar.QueryBegin;
				return IsReverseComplement;
			}
		};
		// define our anchor structure
		struct Anchor {
			string Name;
			unsigned int Length;
			unsigned int Begin;
			unsigned int End;
		};
		// checks to see if this is truly a Mosaik aligned read archive
		static bool CheckMosaikAlignedReadArchive(const char* filename, bool exitOnError);
		// closes the aligned read archive
		void Close(void);
		// loads the next read from the aligned read archive
		bool LoadNextRead(vector<AlignedPosition>& positions, unsigned short& readLength, string& readName);
		// loads the anchors for the aligned read archive
		void LoadAnchors(vector<Anchor>& anchors);
		// returns the number of reads in this file
		mosaik_64u_t GetNumReads(void) const;

	private:
		//
		struct HeaderStatistics {
			off_type AnchorsOffset;
			off_type ReadNamesOffset;
			off_type CurrentReadNameOffset;
			mosaik_64u_t Timestamp;
			mosaik_64u_t NumReads;
			unsigned int NumAnchors;

			HeaderStatistics()
				: AnchorsOffset(0)
				, ReadNamesOffset(0)
				, Timestamp(0)
				, NumReads(0)
				, NumAnchors(0)
			{}
		} mStatistics;
		// checks if the buffer size is large enough to accomodate the requested size
		void CheckBufferSize(const unsigned int requestedBytes);
		// checks if the read buffer size is large enough to accomodate the requested size
		void CheckReadBufferSize(const unsigned int requestedBytes);
		// parses the aligned read archive header
		void ParseHeader(void);
		// denotes whether or not our file has already been closed
		bool mIsStreamOpen;
		// our input file stream
		FILE* mInStream;
		// specifies the buffer used when reading both read and index entries
		unsigned char* mBuffer;
		// specifies the current buffer size
		unsigned int mBufferLen;
		// specifies the buffer used when reading both read and index entries
		char* mReadBuffer;
		// specifies the current buffer size
		unsigned int mReadBufferLen;
		// stores the current read number
		mosaik_64u_t mCurrentReadNum;
	};
}
