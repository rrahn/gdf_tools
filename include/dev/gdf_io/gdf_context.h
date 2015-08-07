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
// Implements the header for the journal sequence file formats.
// ==========================================================================

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_HEADER_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_HEADER_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct SnpCompStrategy2BitEncoding_;
typedef Tag<SnpCompStrategy2BitEncoding_> SnpCompStrategy2BitEncoding;

struct SnpCompStrategyRaw_;
typedef Tag<SnpCompStrategyRaw_> SnpCompStrategyRaw;

typedef TagList<SnpCompStrategy2BitEncoding,
        TagList<SnpCompStrategyRaw> > SnpCompStrategies;

typedef TagSelector<SnpCompStrategies> SnpCompStrategieSelector;


struct GdfIOCompressionContext
{
    CharString                refID;
    CharString                refURI;
    SnpCompStrategieSelector  snpCompStrategy;

    GdfIOCompressionContext
    {
        assign(snpCompStrategy, SnpCompStrategyRaw);
    }

};

// ----------------------------------------------------------------------------
// Class GdfIOContext
// ----------------------------------------------------------------------------

template <typename TNameStore_        = StringSet<CharString>,
          typename TNameStoreCache_   = NameStoreCache<TNameStore_>,
          typename TStorageSpec       = Owner<> >
class GdfIOContext : public GdfIOCompressionContext
{
public:
    typedef TNameStore_         TNameStore;
    typedef TNameStoreCache_    TNameStoreCache;

    typedef typename StorageSwitch<TNameStore, TStorageSpec>::Type      TNameStoreMember;
    typedef typename StorageSwitch<TNameStoreCache, TStorageSpec>::Type TNameStoreCacheMember;

    // Cache for the sequence name lookup.
    TNameStoreMember        _contigNames;
    TNameStoreCacheMember   _contigNamesCache;

    // Cache for the sample name lookup.
    TNameStoreMember        _sampleNames;
    TNameStoreCacheMember   _sampleNamesCache;

    CharString              buffer;

    GdfIOContext() :
        GdfIOCompressionContext(),
        _contigNames(TNameStoreMember()),
        _contigNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(), (TNameStoreCache*)NULL, _contigNames)),
        _sampleNames(TNameStoreMember()),
        _sampleNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(), (TNameStoreCache*)NULL, _sampleNames))
    {}

    GdfIOContext(TNameStore & nameStore_, TNameStoreCache & nameStoreCache_) :
        GdfIOCompressionContext(),
        _contigNames(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(nameStore_)),
        _contigNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(), &nameStoreCache_, _contigNames)),
        _sampleNames(TNameStoreMember()),
        _sampleNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(), (TNameStoreCache*)NULL, _sampleNames))
    {}

    template <typename TOtherStorageSpec>
    GdfIOContext(GdfIOContext<TNameStore, TNameStoreCache, TOtherStorageSpec> & other) :
        GdfIOCompressionContext(),
        _contigNames(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(contigNames(other))),
        _contigNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(), &contigNamesCache(other), _contigNames)),
        _sampleNames(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(sampleNames(other))),
        _sampleNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(), &sampleNamesCache(other), _sampleNames))
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================
}

#endif // INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_HEADER_H_
