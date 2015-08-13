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

#ifndef INCLUDE_DEV_GDF_IO_GDF_CONTEXT_H_
#define INCLUDE_DEV_GDF_IO_GDF_CONTEXT_H_

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

/*!
 * @class GdfIOCompressionContext
 * @headerfile <dev/gdf_io.h>
 * @brief The I/O context to use for GDF I/O.
 */

struct GdfIOCompressionContext
{
    CharString                refID;
    CharString                refUri;
    SnpCompStrategieSelector  snpCompStrategy;

    GdfIOCompressionContext()
    {
        assign(snpCompStrategy, SnpCompStrategyRaw());
    }

};

// ----------------------------------------------------------------------------
// Class GdfIOContext
// ----------------------------------------------------------------------------

/*!
 * @class GdfIOContext
 * @headerfile <dev/gdf_io.h>
 * @brief The I/O context to use for GDF I/O.
 * @implements GdfIOCompressionContext
 *
 * @signature template <typename TNameStore[, typename TNameStoreCache]>
 *            class GdfIOContext : public GdfIOCompressionContext;
 *
 * @tparam TNameStore      The type used to represent the names.
 * @tparam TNameStoreCache The type used to cache the names. Defaults to @link NameStoreCache @endlink &lt;TNameStore&gtl;.
 */

/*!
 * @fn GdfIOContext::GdfIOContext
 * @brief Constructor.
 *
 * @signature GdfIOContext::GdfIOContext();
 * @signature GdfIOContext::GdfIOContext(contigNames, sampleNames);
 *
 * Default constructor or construction with references to contig and sample names.
 */

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

// ----------------------------------------------------------------------------
// Function contigNames()
// ----------------------------------------------------------------------------

/*!
 * @fn GdfIOContext#contigNames
 * @brief Return reference to the contig names from @link GdfIOContext @endlink.
 *
 * @signature TNameStoreRef contigNames(context);
 *
 * @param[in] context The @link GdfIOContext @endlink to query.
 *
 * @return TNameStoreRef A reference to the <tt>TNameStore</tt> of the context.
 */

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore &
contigNames(GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStore &>(context._contigNames);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore const &
contigNames(GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    return _referenceCast<TNameStore &>(context._contigNames);
}

/*!
 * @fn GdfIOContext#contigNamesCache
 * @brief Return reference to contig names cache from @link GdfIOContext @endlink.
 *
 * @signature TNameStoreCacheRef contigNamesCache(context);
 *
 * @param[in] context The @link GdfIOContext @endlink to query.
 *
 * @return TNameStoreCacheRef A reference to the <tt>TNameStoreCache</tt> of the context.
 */

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache &
contigNamesCache(GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStoreCache &>(context._contigNamesCache);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache const &
contigNamesCache(GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    return _referenceCast<TNameStoreCache &>(context._contigNamesCache);
}

/*!
 * @fn GdfIOContext#sampleNames
 * @brief Return reference to the sample names from @link GdfIOContext @endlink.
 *
 * @signature TNameStoreRef sampleNames(context);
 *
 * @param[in] context The @link GdfIOContext @endlink to query.
 *
 * @return TNameStoreRef A reference to the <tt>TNameStore</tt> of the context.
 */

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore &
sampleNames(GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStore &>(context._sampleNames);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore const &
sampleNames(GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    return _referenceCast<TNameStore &>(context._sampleNames);
}

/*!
 * @fn GdfIOContext#sampleNamesCache
 * @brief Return reference to the sample names cache from @link GdfIOContext @endlink.
 *
 * @signature TNameStoreCacheRef sampleNamesChache(context);
 *
 * @param[in] context The @link GdfIOContext @endlink to query.
 *
 * @return TNameStoreCacheRef A reference to the <tt>TNameStoreCache</tt> of the context.
 */

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache &
sampleNamesCache(GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStoreCache &>(context._sampleNamesCache);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache const &
sampleNamesCache(GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    return _referenceCast<TNameStoreCache &>(context._sampleNamesCache);
}
}

#endif // INCLUDE_DEV_GDF_IO_GDF_CONTEXT_H_
