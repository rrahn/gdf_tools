// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Implements the read methods to load journal sequence format.
// ==========================================================================

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_READ_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_READ_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
_parseGdfSamples(GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
                 CharString const & sampleNames)
{
    StringSet<CharString> samples;
    strSplit(samples, sampleNames, EqualsChar<':'>(), false);

    if (empty(samples))
        throw ParseError("Could not determine sample names!");

    for (auto& sampleName : samples)
        nameToId(nameStoreCache(context), sampleName);
}

// ----------------------------------------------------------------------------
// Function readHeader(); Gdf
// ----------------------------------------------------------------------------

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec, typename TForwardIter>
inline void
readHeader(GdfHeader & header,
           GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Gdf /*tag*/)
{
    // Read file information.
    write(buffer, iter, 4);  // Read 4 bytes from input iterator.
    if (*iter != GDF_VERSION_MAJOR)
        throw ParseError("Incompatible mjaor version number!");
    skipOne(iter);
    if (*iter > GDF_VERSION_MINOR)
        throw ParseError("Incompatible minor version number!");
    skipNChars(iter, 2);

    CharString buffer;
    readLine(buffer, iter);

    if (!startsWith(buffer, "!!ref_id="))
    {
        throw ParseError("Expected line starting with \"!!ref_id=\"!");
    }
    else
    {
        context.refID = suffix(buffer, 9);
    }
    clear(buffer);
    readLine(buffer, iter);
    if (!startsWith(buffer, "!!ref_uri="))
    {
        throw ParseError("Expected line starting with \"!!ref_uri=\"!");
    }
    else
    {
        context.refID = suffix(buffer, 10);
    }
    clear(buffer);
    readLine(buffer, iter);
    if (!startsWith(buffer, "!!snp_compression="))
    {
        throw ParseError("Expected line starting with \"!!ref_uri=\"!");
    }
    else
    {
        if (back(buffer) == '2')
            assign(context.snpCompStrategy, SnpCompStrategy2BitEncoding());
        else if (back(buffer) == '1')
            assign(context.snpCompStrategy, SnpCompStrategyRaw());
        else
            throw ParseError("Unsupported snp compression strategy!");
    }

    GdfHeaderRecord record;
    // Read additional information.
    while (!atEnd(iter) && *iter == '#')
    {
        skipOne(iter);
        clear(buffer);

        if (value(iter) == '#')
        {
            // Is header line.
            skipOne(iter);
            clear(record);

            // Read header key.
            readUntil(record.key, iter, OrFunctor<EqualsChar<'='>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gdf> >());

            // Skip '='.
            skipOne(iter);

            // Read header value.
            readLine(record.value, iter);

            // Read sample names;
            if (record.key == "sample_names")
            {
                _parseGdfSamples(context, record.value);
                return;
            }
            appendValue(header, record);
        }
        else
        {
            throw ParseError("Invalid file header!\nExpected two \"##\" at the beginning of this line!\nAre you sure you are reading a gdf file?");
        }
    }
}

template <typename TRecordIt, typename TBufferIt>
inline void
_extractSnp(TRecordIt & recordIt,
            __uint32 & /*deltaPos*/,
            TBufferIt & bufferIt,
            SnpCompStrategyRaw const & /*tag*/)
{
    resize(recordIt->insValue, 1);
    readOne(recordIt->insValue[0], bufferIt);
}

template <typename TRecordIt, typename TBufferIt>
inline void
_extractSnp(TRecordIt & recordIt,
            __uint32 & deltaPos,
            TBufferIt & /*bufferIt*/,
            SnpCompStrategy2BitEncoding const & /*tag*/)
{
    char val = convert<Dna>(deltaPos >> ((sizeof(__uint32) << 3) - 3));
    append(recordIt->insValue, val);
    deltaPos &= (~static_cast<__uint32>(0) >> 3);
    if (isBitSet(deltaPos, 30) || isBitSet(deltaPos, 29))
        throw ParseError("GdfFileIn: Error while reading delta position for snp.");
}

template <typename TBuffer, typename TSnpCompression>
inline void
_extractDataBlock(GdfRecord & record,
                  TBuffer const & buffer,
                  TSnpCompression const & /*tag*/)
{
    auto bufferIt = begin(buffer);
    auto recordIt = begin(record);

    __uint32 lastRefPos = 0;
    // Write the reference offset of the current block.
    readRawPod(lastRefPos, bufferIt);
    for (; recordIt != end(record); ++recordIt)
    {
        clear(recordIt->insValue);
        recordIt->contigPos = lastRefPos;
        __uint32 deltaPos = 0;
        readRawPod(deltaPos, bufferIt);
        // Read a snp.
        if (isBitSet(deltaPos, BitsPerValue<__uint32>::VALUE - 1))
        {
            clearBit(deltaPos, BitsPerValue<__uint32>::VALUE - 1);
            _extractSnp(recordIt, deltaPos, bufferIt, TSnpCompression());
            recordIt->deltaType = DELTA_TYPE_SNP;
        }
        else  // Read Indel or SV.
        {
            __uint32 valSize;
            readRawPod(valSize, bufferIt);
            if (isBitSet(valSize, 31))  // Extract deletion.
            {
                clearBit(valSize, 31);
                recordIt->delValue = valSize;
                recordIt->deltaType = DELTA_TYPE_DEL;
            }
            else  // This is an insertion or a SV.
            {
                readRawPod(valSize, bufferIt);
                if (isBitSet(valSize, 30)) // This is a SV.
                {
                    clearBit(valSize, 30);
                    write(recordIt->insValue, bufferIt, valSize);
                    readRawPod(recordIt->delValue, bufferIt);
                    recordIt->deltaType = DELTA_TYPE_SV;
                }
                else  // Read insertion.
                {
                    write(recordIt->insValue, bufferIt, valSize);
                    recordIt->deltaType = DELTA_TYPE_INS;
                }
            }
        }
        recordIt->contigPos = lastRefPos + deltaPos;
        lastRefPos = recordIt->contigPos;
    }

    // Now we read the coverage!

    __uint32 hostSize;
    readRawPod(hostSize, bufferIt);
    recordIt = begin(record);
    for (; recordIt != end(record); ++recordIt)
        write(recordIt->coverage, bufferIt, hostSize);
}


template <typename TBuffer>
inline void
_extractDataBlock(GdfRecord const & /*record*/,
                  TBuffer const & /*buffer*/,
                  TagSelector<> const & /*format*/)
{
    throw ParseError("GdfFileIn: Could not determine correct snp compression strategy!");
}

template <typename TBuffer, typename TTagList>
inline void
_extractDataBlock(GdfRecord & record,
                  TBuffer & buffer,
                  TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        _extractDataBlock(record, buffer, TFormat());
    else
        _extractDataBlock(record, buffer, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// ----------------------------------------------------------------------------
// Function readRecord()                                            [GdfRecord]
// ----------------------------------------------------------------------------
// Read record, updating list of known sequences if new one occurs.

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(GdfRecord & record,
           GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Gdf const & /*tag*/)
{
    __uint32 recordNum = 0;
    __uint32 blockSize = 0;
    readRawPod(recordNum, iter);
    resize(record, recordNum, Exact());
    clear(record);
    readRawPod(blockSize, iter);
    clear(context.buffer);
    write(context.buffer, iter, blockSize);
    _extractDataBlock(record, context.buffer, context.snpCompStrategy);
}

}

#endif // INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_READ_H_
