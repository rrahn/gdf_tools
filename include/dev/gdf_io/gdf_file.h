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
// Formatted File definition for GDF file.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GDF_IO_GDF_FILE_H_
#define INCLUDE_SEQAN_GDF_IO_GDF_FILE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Typedef GdfFileIn
// ----------------------------------------------------------------------------

/*!
 * @class GdfFileIn
 * @signature typedef FormattedFile<Gdf, Input> GdfFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/gdf_io.h>
 * @brief Class for reading GDF files.
 *
 * @see GffRecord
 */

typedef FormattedFile<Gdf, Input>   GdfFileIn;

// ----------------------------------------------------------------------------
// Typedef GffFileOut
// ----------------------------------------------------------------------------

/*!
 * @class GdfFileOut
 * @signature typedef FormattedFile<Gdf, Output> GdfFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/gdf_io.h>
 * @brief Class for writing GDF files.
 *
 * @see GffRecord
 */

typedef FormattedFile<Gdf, Output>  GdfFileOut;

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Gdf, T>
{
    static unsigned char const VALUE[4];
};

template <typename T>
unsigned char const MagicHeader<Gdf, T>::VALUE[4] =
{
    '#', 'G', 'D', 'F'  // GDF's magic header
};

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Gdf, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Gdf, T>::VALUE[1] =
{
    ".gdf"     // default output extension
};


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TSpec, typename TDirection, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Gdf, TDirection, TSpec>, TStorageSpec>
{
    typedef StringSet<CharString>                                   TNameStore;
    typedef NameStoreCache<TNameStore>                              TNameStoreCache;
    typedef GdfIOContext<TNameStore, TNameStoreCache, TStorageSpec> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<Gdf, TDirection, TSpec> >
{
    typedef TagSelector<
                TagList<Gdf>
            > Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readHeader(); GdfHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readHeader(GdfHeader & header, FormattedFile<Gdf, Input, TSpec> & file)
{
    readHeader(header, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

// convient GdfFile variant
template <typename TSpec>
inline void
readRecord(GdfRecord & record, FormattedFile<Gdf, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function writeHeader(); GdfHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeHeader(FormattedFile<Gdf, Output, TSpec> & file, GdfHeader & header)
{
    writeHeader(file.iter, header, context(file), file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(); GdfRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(FormattedFile<Gdf, Output, TSpec> & file, GdfRecord & record)
{
    writeRecord(file.iter, record, context(file), file.format);
}

}

#endif // INCLUDE_SEQAN_GDF_IO_GDF_FILE_H_
