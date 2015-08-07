// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements the write methods to write out the journal strings as
// binary version.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GDF_IO_GDF_WRITE_H_
#define INCLUDE_SEQAN_GDF_IO_GDF_WRITE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

static const char DIFF_FILE_SEPARATOR = '\t';

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function writeHeader(); Gdf
// ----------------------------------------------------------------------------

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
writeHeader(TTarget & target,
            GdfHeader const & header,
            GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
            Gdf /*tag*/)
{
    // Write file information.
    write(target, &MagicHeader<Gdf>::VALUE[0], 4);
    writeValue(GDF_VERSION_MAJOR);
    writeValue(GDF_VERSION_MINOR);

    // Write reference file information.
    if (empty(context.refID))
        throw (ParseError("Reference ID not set!"));
    if (empty(context.refURI))
        throw (ParseError("Reference URI not set!"));

    write(target, "!!ref_id=");
    write(target, context.refID);
    writeValue(target, '\n');
    write(target, "!!ref_uri=");
    write(target, context.refURI);
    writeValue(target, '\n');

    // Write compression information.
    write(target, "!!snp_compression=");
    if (isEqual(context.snpCompStrategy, SnpCompStrategy2BitEncoding))
        writeValue(target, '2');
    else if (isEqual(context.snpCompStrategy, SnpCompStrategy2BitEncoding))
        writeValue(target, '1');
    else
        throw (ParseError("Unknown snp compression strategy!"));
    writeValue(target, '\n');

    // Write additional information.
    for (auto& headerRecord : header)
    {
        write(target, "##");
        write(target, headerRecord.key);
        writeValue(target, '=');
        write(target, headerRecord.value);
        writeValue(target, '\n');
    }

    // Write sample names.
    write(target, "##sample_names=")
    for (auto& name : sampleNames(context))
    {
        write(target, name);
        writeValue(target, ':');
    }
    writeValue(target, '\n');
}

// ----------------------------------------------------------------------------
// Function _writeBitVector()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TSpec>
inline void
_writeBitVector(TTarget & target, String<bool, Packed<TSpec> > const & bitVec)
{
    typedef String<bool, Packed<TSpec> > const TConstBitVec;
    typedef typename Host<TConstBitVec>::Type THost;
    typedef typename Value<THost>::Type THostValue;

    THostValue* hostVal = &(host(bitVec)[0]);
    write(target, hostVal, length(host(bitVec)));
//    for (; it != itEnd; ++it)
//        appendRawPod(target, reinterpret_cast<const char *>(&it->i), sizeof(TBitVector));
}

// ----------------------------------------------------------------------------
// Function _writeSnp()
// ----------------------------------------------------------------------------

// Writes SNP in separat value after delta pos.
template <typename TTarget, typename TPosition, typename TAlphabet>
inline void
_writeSnp(TTarget & buffer,
          TPosition deltaPos,
          TAlphabet snp,
          SnpCompStrategyRaw /*tag*/)
{
    SEQAN_ASSERT_NOT(isBitSet(deltaPos, BitsPerValue<__uint32>::VALUE -1));  // The last bit should not be set.

    setBit(deltaPos, BitsPerValue<__uint32>::VALUE -1);  // Set bit to indicate SNP
//    register __uint32 offset = sizeof(__uint32) + sizeof(TAlphabet);
//    char* refBuffer = reinterpret_cast<char*>(&deltaPos);
//    resize(blockBuffer, length(blockBuffer) + offset);
//
//    if (SystemByteOrder::IS_LITTLE_ENDIAN())
//        arrayMoveForwardReverse(refBuffer, refBuffer + sizeof(__uint32), end(blockBuffer, Standard()) - offset);
//    else
//        arrayMoveForward(refBuffer, refBuffer + sizeof(__uint32), end(blockBuffer, Standard()) - offset);
    appendRawPod(buffer, deltaPos);
    writeValue(buffer, snp);
//    char* snpBuffer = reinterpret_cast<char*>(&snp);
//    arrayMoveForward(snpBuffer, snpBuffer + sizeof(TAlphabet), end(blockBuffer, Standard()) - sizeof(TAlphabet));
}

// ----------------------------------------------------------------------------
// Function _writeSnp()                                                   [Dna]
// ----------------------------------------------------------------------------

// Encodes SNP in delta position assuming, that 28 bits are sufficient to store the delta to the previous variant.
template <typename TTarget, typename TPosition, typename TAlphabet>
inline void
_writeSnp(TTarget & buffer,
          TPosition deltaPos,
          TAlphabet snp,
          SnpCompStrategy2BitEncoding /*tag*/)
{
    SEQAN_ASSERT_NOT(isBitSet(deltaPos, BitsPerValue<__uint32>::VALUE -1));  // Bit at index 31 not set.
    SEQAN_ASSERT_NOT(isBitSet(deltaPos, BitsPerValue<__uint32>::VALUE -2));  // Bit at index 30 not set.
    SEQAN_ASSERT_NOT(isBitSet(deltaPos, BitsPerValue<__uint32>::VALUE -3));  // Bit at index 29 not set.

    setBit(deltaPos, BitsPerValue<__uint32>::VALUE -1);
    deltaPos |= static_cast<__uint32>(convert<Dna>(snp)) << ((sizeof(__uint32) << 3) - 3);
//    resize(blockBuffer, length(blockBuffer) + sizeof(__uint32));
//    char* refBuffer = reinterpret_cast<char*>(&deltaPos);
//    if (SystemByteOrder::IS_LITTLE_ENDIAN())
//        arrayMoveForwardReverse(refBuffer, refBuffer + sizeof(__uint32), end(blockBuffer, Standard()) - sizeof(__uint32));
//    else
//        arrayMoveForward(refBuffer, refBuffer + sizeof(__uint32), end(blockBuffer, Standard()) - sizeof(__uint32));
    appendRawPod(buffer, deltaPos);
}

// ----------------------------------------------------------------------------
// Function _writeDataBlock()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TRecordIter, typename TSnpCompressionStrategy>
inline int _writeDataBlock(TTarget & target,
                           TRecordIter & itBegin,
                           TRecordIter & itEnd,
                           TSnpCompressionStrategy const & /*deltaMap*/)
{
    typedef DeltaMap<TValue, TAlphabet> TDeltaMap;

    CharString blockBuffer;
    typename Iterator<CharString>::Type blockBuffIt;
    __uint32 lastRefPos = itBegin->contigPos;
    // Write the reference offset of the current block.
    appendRawPod(target, lastRefPos);

    TRecordIter it = itBegin;
    while(it != itEnd)
    {
        // Store the reference position.
        __uint32 deltaPos = it->contigPos - lastRefPos;

        // Write SNP Data.
        if (it->deltaType == DeltaType::DELTA_TYPE_SNP)
        {
//            resize(blockBuffer, length(blockBuffer) + sizeof(__uint32));
//            blockBuffIt = end(blockBuffer) - sizeof(__uint32);
            _writeSnp(target, deltaPos, it->insValue[0], TSnpCompressionStrategy());
//            setBit(deltaPos, BitsPerValue<__uint32>::VALUE -1);
//            deltaPos |= static_cast<__uint32>(deltaSnp(deltaMap, deltaPosition(deltaInfo))) << (sizeof(__uint32) * 8 - 3);
//            char* refBuffer = reinterpret_cast<char*>(&deltaPos);
//
//            if (SystemByteOrder::IS_LITTLE_ENDIAN())
//                arrayMoveForwardReverse(refBuffer, refBuffer + sizeof(__uint32), blockBuffIt);
//            else
//                arrayMoveForward(refBuffer, refBuffer + sizeof(__uint32), blockBuffIt);
        }
        else  // Write Indel
        {
//            char * refBuffer = reinterpret_cast<char *>(&deltaPos);
//            resize(blockBuffer, length(blockBuffer) + sizeof(__uint32));
//            blockBuffIt = end(blockBuffer) - sizeof(__uint32);
//            if (SystemByteOrder::IS_LITTLE_ENDIAN())
//                arrayMoveForwardReverse(refBuffer, refBuffer + sizeof(__uint32), blockBuffIt);
//            else
//                arrayMoveForward(refBuffer, refBuffer + sizeof(__uint32), blockBuffIt);
            appendRawPod(target, deltaPos);

            if (it->deltaType == DeltaType::DELTA_TYPE_DEL)  // Write deletion.
            {
                __uint32 del = static_cast<__uint32>(it->delValue);
                setBit(del, 31);
                appendRawPod(target, del);
//                const char * delBuffer = reinterpret_cast<const char *>(&del);
//                resize(blockBuffer, length(blockBuffer) + sizeof(del));
//                blockBuffIt = end(blockBuffer) - sizeof(del);
//                if (SystemByteOrder::IS_LITTLE_ENDIAN())
//                    arrayMoveForwardReverse(delBuffer, delBuffer + sizeof(del), blockBuffIt);
//                else
//                    arrayMoveForward(delBuffer, delBuffer + sizeof(del), blockBuffIt);
            }
            else
            {
//                typedef typename DeltaValue<TDeltaMap const, DeltaType::DELTA_TYPE_INS>::Type TIns;

                // Handle Indel.
                if (it->deltaType == DeltaType::DELTA_TYPE_INDEL)
                {
//                    TIns ins = deltaIndel(itDelta).i2;
                    __uint32 insLength = length(it->insValue);
                    setBit(insLength, 30);
                    __uint32 del = static_cast<__uint32>(it->delValue);

//                    const char * delBuffer = reinterpret_cast<const char *>(&del);
//                    const char * insBuffer = reinterpret_cast<const char *>(&insLength);
//                    resize(blockBuffer, length(blockBuffer) + sizeof(insLength) + sizeof(del) + length(ins));
//                    blockBuffIt = end(blockBuffer) - sizeof(insLength) - sizeof(del) - length(ins);
//
//                    if (SystemByteOrder::IS_LITTLE_ENDIAN())
//                        arrayMoveForwardReverse(insBuffer, insBuffer + sizeof(insLength), blockBuffIt);
//                    else
//                        arrayMoveForward(insBuffer, insBuffer + sizeof(insLength), blockBuffIt);
                    // Write size of insertion.
                    appendRawPod(target, insLength);
                    // Write inserted characters.
                    write(target, it->insValue);
//                    blockBuffIt = end(blockBuffer) - length(ins) - sizeof(del);
//                    arrayMoveForward(begin(ins, Standard()), end(ins, Standard()), blockBuffIt);
                    // Write size of deletion
                    appendRawPod(target, del);
//                    blockBuffIt = end(blockBuffer) - sizeof(del);
//                    if (SystemByteOrder::IS_LITTLE_ENDIAN())
//                        arrayMoveForwardReverse(delBuffer, delBuffer + sizeof(del), blockBuffIt);
//                    else
//                        arrayMoveForward(delBuffer, delBuffer + sizeof(del), blockBuffIt);
                }
                else  // Handle Insertion.
                {
                    SEQAN_ASSERT(it->deltaType == DeltaType::DELTA_TYPE_INS);

                    __uint32 insLength = length(it->insValue);
                    SEQAN_ASSERT_NOT(isBitSet(insLength, BitsPerValue<__uint32>::VALUE - 2));
                    appendRawPod(target, insLength);
                    write(target, it->insValue);
//
//                    const char * insBuffer = reinterpret_cast<const char *>(&insLength);
//                    resize(blockBuffer, length(blockBuffer) + sizeof(insLength) + length(ins));
//                    blockBuffIt = end(blockBuffer) - sizeof(insLength) - length(ins);
//                    if (SystemByteOrder::IS_LITTLE_ENDIAN())
//                        arrayMoveForwardReverse(insBuffer, insBuffer + sizeof(insLength), blockBuffIt);
//                    else
//                        arrayMoveForward(insBuffer, insBuffer + sizeof(insLength), blockBuffIt);
//                    blockBuffIt = end(blockBuffer) - length(ins);
//                    arrayMoveForward(begin(ins, Standard()), end(ins, Standard()), blockBuffIt);
                }
            }
        }
        lastRefPos = it->contigPos;
        ++it;
    }

    // Write the block Data.
//    unsigned blockLength = length(blockBuffer);
//    streamWriteBlock(stream, reinterpret_cast<char*>(&blockLength), sizeof(blockLength));
//    streamWriteBlock(stream, &blockBuffer[0], blockLength);

    auto it = itBegin;
    typedef typename std::remove_reference<decltype(host(it->coverage)[0])>::type THostValue;
    __uint32 hostSize = length(host(it->coverage)) * sizeof(THostValue());
    appendRawPod(target, hostSize);
    // Write the coverage of the block
    for (; it != itEnd; ++it)
        _writeBitVector(target, it->coverage);

    return 0;
}

// Write the io context.
//template <typename TStream, typename TDeltaStore, typename TDeltaCoverageStore, typename TSize>
//inline int _writeGdfData(TStream & stream,
//                          DeltaMap<TDeltaStore, TDeltaCoverageStore> const & deltaMap,
//                          TSize const & blockSize)
//{
//    typedef DeltaMap<TDeltaStore, TDeltaCoverageStore> TDeltaMap;
//    typedef typename Iterator<TDeltaMap const, Standard>::Type TIterator;
//
//    unsigned maxNumOfNodes = length(deltaMap);
//    unsigned numOfBlocks = std::ceil(static_cast<double>(maxNumOfNodes)/static_cast<double>(blockSize));
//
//    // Write the block containing the number of blocks to read.
//    streamWriteBlock(stream, reinterpret_cast<char*>(&numOfBlocks), sizeof(numOfBlocks));
//
//    for (unsigned i = 0; i < numOfBlocks; ++i)  // For each block:
//    {
//        TIterator it = begin(deltaMap, Standard()) + (blockSize * i);
//        TIterator itEnd = _min(end(deltaMap, Standard()), it + blockSize);
//        _writeDataBlock(stream, it, itEnd, deltaMap);
//    }
//    return 0;
//}

//template <typename TStream, typename TValue, typename TAlphabet, typename TConfig>
//inline int
//write(TStream & stream,
//      DeltaMap<TValue, TAlphabet> const & deltaMap,
//      GdfHeader<TConfig> const & gdfHeader,
//      Gdf const & /*tag*/)
//{
//    int res = _writeGdfHeader(stream, gdfHeader);
//    if (res != 0)
//        return res;
//    return _writeGdfData(stream, deltaMap, gdfHeader._fileInfos._blockSize);
//}

template <typename TBuffer>
inline void
_writeDataBlock(TBuffer & /*buffer*/,
                GdfRecord const & /*record*/,
                TagSelector<> const & /*format*/)
{
    SEQAN_FAIL("GdfFileOut: File format not specified.");
}

template <typename TBuffer, typename TTagList>
inline void
_writeDataBlock(TBuffer & buffer,
                GdfRecord const & record,
                TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        _writeDataBlock(buffer, begin(record.data), end(record.data), TFormat());
    else
        _writeDataBlock(buffer, record, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

template <typename TTarget, typename TContext>
inline void
writeRecord(TTarget & target,
            GdfRecord const & record,
            TContext & ioContext,
            Gdf /*tag*/)
{
    clear(ioContext.buffer);
    // Write the number of entries for this block.
    appendRawPod(ioContext.buffer, static_cast<__int32>(length(record)));
    _writeDataBlock(ioContext.buffer, record, ioContext.snpCompStrategy);
    appendRawPod(target, static_cast<__int32>(length(ioContext.buffer)));
    write(target, &(ioContext.buffer[0]), length(ioContext.buffer));  // Write buffer into stream.
}

///*!
// * @fn journalSet#write
// * @deprecated
// */
//template <typename TStream, typename TJournalSequence>
//inline int
//write(TStream & stream,
//      StringSet<TJournalSequence, Owner<JournaledSet> > const & journalSet,
//      GdfHeader const & jseqHeader,
//      Gdf const & /*tag*/)
//{
//    typedef typename Value<TJournalSequence>::Type TAlphabet;
//    typedef typename Position<TJournalSequence>::Type TPosition;
//
//    if (empty(host(journalSet)))
//        return -1;
//    DeltaMap<TPosition, TAlphabet> deltaMap;
//    adaptTo(deltaMap, journalSet);
//    write(stream, deltaMap, jseqHeader, Gdf());
//    return 0;
//}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_GDF_IO_GDF_WRITE_H_
