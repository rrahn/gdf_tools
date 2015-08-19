// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Implements strategy to adapt a vcf file to journaled data.
// ==========================================================================

#ifndef EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_ADAPT_VCF_H_
#define EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_ADAPT_VCF_H_

#include "gdf_tools.h"
#include <seqan/vcf_io.h>

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

static const char * VCF_RECORD_FORMAT_GENOTYPE = "GT";
static const char   VCF_RECORD_GENOTYPE_SEPARATOR_PHASED = '|';
static const char   VCF_RECORD_GENOTYPE_SEPARATOR_UNPHASED = '/';
static const char   VCF_RECORD_ALT_SYMBOLIC_ALLEL_BEGIN = '<';
static const char   VCF_RECORD_ALT_SYMBOLIC_ALLEL_END = '>';
static const char   VCF_RECORD_ALT_BREAKEND_REPL_DEFAULT = '[';
static const char   VCF_RECORD_ALT_BREAKEND_REPL_REVCOMP = ']';


struct BasicDeltaValue
{
    DeltaType       deltaType;
    CharString      insBases;
    unsigned        delBases;
    unsigned        pos;
    String<unsigned> coverage;
};


template <typename T>
struct VcfRecordTranslator{};

template <typename TSequence, typename TConfig, typename TSpec>
struct VcfRecordTranslator<JournaledStringTree<TSequence, TConfig, TSpec> >
{
    typedef JournaledStringTree<TSequence, TConfig, TSpec> TJst;

    TJst & _jst;

    VcfRecordTranslator(TJst & jst) : _jst(jst)
    {}

    inline void operator()(BasicDeltaValue const & basicDelta, bool force = false)
    {
        if (!force && empty(basicDelta.coverage))
            return;

        switch(basicDelta.deltaType)
        {
            case DELTA_TYPE_SNP:
                insert(_jst, basicDelta.pos, basicDelta.insBases[0], basicDelta.coverage, DeltaTypeSnp());
                break;
            case DELTA_TYPE_INS:
                insert(_jst, basicDelta.pos, basicDelta.insBases, basicDelta.coverage, DeltaTypeIns());
                break;
            case DELTA_TYPE_DEL:
                insert(_jst, basicDelta.pos, basicDelta.delBases, basicDelta.coverage, DeltaTypeDel());
                break;
            default:
                throw(RuntimeError("Unknown delta event in vcf file."));
        }
    }
};

template <typename TPos, typename TSize, typename TInsBuffer>
struct VariantInfo
{
    TPos refPos;
    TSize variantSize;
    DeltaType variantType;
    TInsBuffer seqBuffer;

    VariantInfo() : refPos(0), variantSize(0)
    {}
};


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _resolveConflicts()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TConfig, typename TSpec>
inline void _resolveConflicts(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    typedef JournaledStringTree<TSequence, TConfig, TSpec> TJst;
    typedef typename Member<TJst, JstDeltaMapMember>::Type TDeltaMap;
    typedef typename TDeltaMap::TDeltaEntries TDeltaEntries;
    typedef typename Value<TDeltaEntries>::Type TDeltaEntry;
    typedef typename DeltaRecord<TDeltaEntry>::Type TDeltaRecord;
    typedef typename Iterator<TDeltaMap, Standard>::Type TStoreIter;
    typedef typename DeltaCoverage<TDeltaMap>::Type TBitVector;
    typedef typename Position<TDeltaMap>::Type TPosition;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeIns>::Type TIns;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSV>::Type TSV;

    TDeltaMap& deltaMap = impl::member(jst, JstDeltaMapMember());
    TStoreIter itBegin = begin(deltaMap, Standard());
    TStoreIter it = itBegin;
    TStoreIter itEnd = end(deltaMap, Standard());
    for (;it != itEnd; ++it)
    {
        if (deltaType(it) == DELTA_TYPE_DEL)
        {  // Resolve all variants that intersect with a deleted region.
            TPosition endPoint = deltaValue(it, DeltaTypeDel()) + getDeltaPosition(*it);
            // Move to the begin of the affected range.
            TStoreIter itLocal = it;
            while (itLocal != begin(deltaMap, Standard()) && getDeltaPosition(*(--itLocal)) == getDeltaPosition(*it))
            {}

            if (getDeltaPosition(*itLocal) != getDeltaPosition(*it))  // Move one up to begin of range.
                ++itLocal;

            // Resolve the conflicts.
            while (itLocal != end(deltaMap, Standard()) && getDeltaPosition(*itLocal) < endPoint)
            {
                if (itLocal != it)
                    transform(getDeltaCoverage(*itLocal), getDeltaCoverage(*itLocal), getDeltaCoverage(*it),
                              FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
                ++itLocal;
            }
        }
        if (deltaType(it) == DELTA_TYPE_INS)
        {  // Resolve all insertions that occur immediately before an replacement or deletion.
            TStoreIter itLocal = it + 1;
            while (getDeltaPosition(*itLocal) == getDeltaPosition(*it))
            {
//                TMappedDelta deltaInfoInner = mappedDelta(deltaMap, itLocal - itBegin);
                if (deltaType(itLocal) == DELTA_TYPE_DEL || deltaType(itLocal) == DELTA_TYPE_SNP)
                {
                    TBitVector tmpVec;
                    transform(tmpVec, getDeltaCoverage(*it), getDeltaCoverage(*itLocal), FunctorBitwiseAnd());
                    if (!testAllZeros(tmpVec))  // At least one sequence contains an ambiguous variant.
                    {  // TODO(rrahn): Check if all variants are reset.
                        // Update the coverage of the nodes.
                        transform(getDeltaCoverage(*it), getDeltaCoverage(*it), tmpVec,
                                  FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
                        transform(getDeltaCoverage(*itLocal), getDeltaCoverage(*itLocal), tmpVec,
                                  FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());

                        if (deltaType(itLocal) == DELTA_TYPE_DEL)
                        {
                            appendValue(getDeltaStore(deltaMap._deltaStore, DeltaTypeSV()),
                                                      TSV(deltaValue(itLocal, DeltaTypeDel()),
                                                          deltaValue(it, DeltaTypeIns())));
                        }
                        else
                        {
                            TIns tmp = deltaValue(it, DeltaTypeIns());
                            append(tmp, deltaValue(itLocal, DeltaTypeSnp()));
                            appendValue(getDeltaStore(deltaMap._deltaStore, DeltaTypeSV()), TSV(1, tmp));
                        }
                        // Insert the new coverage and variant info into the variant store.
                        TPosition currPos = it - begin(deltaMap, Standard());
                        insertValue(deltaMap._entries, currPos,
                                    TDeltaEntry(getDeltaPosition(*it), TDeltaRecord(DELTA_TYPE_SV,
                                                length(getDeltaStore(deltaMap._deltaStore, DeltaTypeSV())) - 1),
                                                tmpVec));
//                        insertValue(deltaMap._deltaCoverageStore._coverageData, currPos, tmpVec);
//                        insertValue(deltaMap._deltaStore._varDataMap, currPos,
//                                    (length(deltaMap._deltaStore._indelData) -1) | DELTA_TYPE_SV);

                        // Insert the ref position and synchronize the iterator.
//                        insertValue(deltaMap._keys, currPos, *it);
                        it = begin(deltaMap, Standard()) + currPos + 1;
                        itEnd = end(deltaMap, Standard());  // Update end pointer.
                        itLocal = it +1;
                    }
                }
                ++itLocal;
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function _extracVariants()
// ----------------------------------------------------------------------------

//template <typename TVariantInfo, typename TVcfRecord>
//inline int _extractVariants(String<TVariantInfo> & varString,
//                            TVcfRecord & vcfRecord)
//{
//
//    StringSet<CharString> altVarSet;
//    splitString(altVarSet, vcfRecord.alt, ',');
//
//    resize(varString, length(altVarSet), Exact());
//
//    for (unsigned i = 0; i < length(varString); ++i)  // Iterate over the variants.
//    {
//        CharString alt = altVarSet[i];
//
//        int leftLcp = lcpLength(vcfRecord.ref, alt);
//        CharString altR(suffix(alt, leftLcp));  // Only take from left lcp different value.
//        CharString refR(suffix(vcfRecord.ref, leftLcp));  // Only take from right lcp different value.
//        reverse(altR);
//        reverse(refR);
//        int rightLcp = lcpLength(altR, refR);
//
//        // Extract the corresponding changes.
//        Segment<CharString, InfixSegment> refInf = infix(vcfRecord.ref, leftLcp, length(vcfRecord.ref) - rightLcp);
//        Segment<CharString, InfixSegment> altInf = infix(alt, leftLcp, length(alt) - rightLcp);
//
//        TVariantInfo & info = varString[i];
//        info.refPos = vcfRecord.beginPos + leftLcp;
//
//        // Variant is a SNP or a single base deletion.
//        if (length(refInf) == length(altInf))
//        {
//            if (length(alt) == 1 && alt[0] == '.')  // Single base deletion
//            {
//    //            std::cerr << "DEL ref: " << refInf << " alt: " << altInf << std::endl;
//                info.variantType = DELTA_TYPE_DEL;
//                info.variantSize = length(refInf);
//            }
//            else
//            {
//    //            std::cerr << "SNP ref: " << refInf << " alt: " << altInf << std::endl;
//                info.variantType = DELTA_TYPE_SNP;
//                info.variantSize = length(refInf);
//                info.seqBuffer = altInf;
//            }
//        }
//        // Variant is a deletion.
//        else if (length(altInf) == 0)
//        {
//    //        std::cerr << "DEL ref: " << refInf << " alt: " << altInf << std::endl;
//            info.variantType = DELTA_TYPE_DEL;
//            info.variantSize = length(refInf);
//        }
//        // Variant is an Insertion
//        else if (length(refInf) == 0)
//        {
//    //        std::cerr << "INS ref: " << refInf << " alt: " << altInf << std::endl;
//            info.variantType = DELTA_TYPE_INS;
//            info.variantSize = length(altInf);
//            info.seqBuffer = altInf;
//        }
//
//        // TODO(rrahn): Use structural variants.
//    }
//    return 0;
//}

// ----------------------------------------------------------------------------
// Function _extractHaplotype()
// ----------------------------------------------------------------------------

//template <typename TSource>
//inline unsigned
//_extractHaplotype(TSource const & source)
//{
//    // Fast switch for diploides.
//    switch(source)
//    {
//        case '0': return 0;
//        case '1': return 1;
//        case '2': return 2;
//        default:
//            unsigned res;
//            std::stringstream buffer;
//            buffer << source;
//            buffer >> res;
//            return res;
//    }
//}
//
//inline unsigned
//numOfHaploytpesPerSequence(VcfHeader const & header)
//{
//    for (auto& headerRecord : header)
//    {
//        if (headerRecord.key == VCF_HEADER_PLOIDY_KEY)
//            return lexicalCast<unsigned>(headerRecord.value);
//    }
//    return 1;  // Assume monoploid if no ploidy information is set.
//}

inline bool
isGenotypePresent(VcfRecord const & record)
{
    return startsWith(record.format, VCF_RECORD_FORMAT_GENOTYPE);
}

inline String<String<unsigned> >
extractGenotypeInfos(VcfRecord const & record, StringSet<CharString> const & altSet)
{
    typedef OrFunctor<EqualsChar<VCF_RECORD_GENOTYPE_SEPARATOR_PHASED>, EqualsChar<VCF_RECORD_GENOTYPE_SEPARATOR_UNPHASED> > TPhasingSeparator;
    typedef AssertFunctor<OrFunctor<IsDigit, EqualsChar<'.'> >, ParseError> TAssertFunctor;

    String<String<unsigned> > altCoverageSet;
    resize(altCoverageSet, length(altSet), Exact());

    // The allel information (GT field in VCF) must be present and always comes at the first part.
    if (isGenotypePresent(record))
    {
        // Parse genotype field per sample.
        for (auto indIt = begin(record.genotypeInfos, Standard()); indIt != end(record.genotypeInfos, Standard()); ++indIt)
        {
            StringSet<CharString> genotypeField;
            strSplit(genotypeField, *indIt, EqualsChar<':'>());

            if (empty(front(genotypeField)))
                continue;  // continue with next individual.

            auto allelIt = begin(front(genotypeField), Standard());
            CharString buffer;
            String<int> altIds;
            do
            {
                clear(buffer);
                readUntil(buffer, allelIt, TPhasingSeparator(), TAssertFunctor());
                skipOne(allelIt);

                if (*allelIt == VCF_RECORD_GENOTYPE_SEPARATOR_UNPHASED)  // We don't support unphased genotypes
                {
                    clear(altIds);
                    break;
                }

                SEQAN_ASSERT_EQ(length(buffer), 1u);

                if (buffer == ".")  // ignore unknown fields.
                    appendValue(altIds, -1);

                // Extract the haplotype specific position.
                appendValue(altIds, lexicalCast<int>(buffer) - 1);
            } while (allelIt != end(front(genotypeField)));

            auto ploidy = length(altIds);
//            decltype(ploidy) haplotypeNum = 0;

            for (decltype(ploidy) htNum = 0; htNum < length(altIds); ++htNum)
            {
                if (altIds[htNum] >= 0)
                    appendValue(altCoverageSet[altIds[htNum]], ploidy * (indIt - begin(record.genotypeInfos, Standard())) + htNum);
            }
        }
    }
    return altCoverageSet;
}

template <typename TString>
inline bool
isSymbolicAllel(TString const & alt)
{
    if (front(alt) != VCF_RECORD_ALT_SYMBOLIC_ALLEL_BEGIN)
        return false;

    if (back(alt) == VCF_RECORD_ALT_SYMBOLIC_ALLEL_END)
        return true;

    CharString message = "Unknown alternative allel description: ";
    append(message, alt);

    throw(ParseError(toCString(message)));
    return false;
}

template <typename TString>
inline bool
isBreakendReplacement(TString const & alt)
{
    if (front(alt) == VCF_RECORD_ALT_BREAKEND_REPL_DEFAULT || front(alt) == VCF_RECORD_ALT_BREAKEND_REPL_REVCOMP)
        return true;
    if (back(alt) == VCF_RECORD_ALT_BREAKEND_REPL_DEFAULT || back(alt) == VCF_RECORD_ALT_BREAKEND_REPL_REVCOMP)
        return true;
    return false;
}

inline BasicDeltaValue
extractBasicDeltaValue(VcfRecord const & record, CharString const & alt)
{
    typedef ModifiedString<CharString const, ModReverse> TModReverse;

    // Sanity check.
    SEQAN_ASSERT(!isSymbolicAllel(alt));
    SEQAN_ASSERT(!isBreakendReplacement(alt));

    if (empty(record.ref) || empty(alt))
        throw(IOError("Invalid REF or ALT value in VCF format."));

    TModReverse revRef(record.ref);
    TModReverse revAlt(alt);

    auto lcp = lcpLength(record.ref, alt);
    auto lcs = lcpLength(revRef, revAlt);
    BasicDeltaValue val;
    val.pos = record.beginPos;
    if (length(record.ref) == length(alt))  // Simple replacement.
    {
        SEQAN_ASSERT_LT(lcp + lcs, length(record.ref));  // The lcp + lcs must be less than the actual length of the ref.
        SEQAN_ASSERT_LT(lcp, length(alt) - lcs);
        if (length(record.ref) == 1)
            val.deltaType = DELTA_TYPE_SNP;
        else
            val.deltaType = DELTA_TYPE_SV;
            val.insBases = infix(alt, lcp, length(alt) - lcs);
        val.delBases = length(val.insBases);
        --val.pos;
    }
    else  // Extract small InDel
    {
        bool refBaseAfter = (record.beginPos == 1);  // REF base is the one after the event.

        if (length(record.ref) < length(alt))  // Simple insertion.
        {
            if (refBaseAfter)
            {
                val.insBases = prefix(alt, length(alt) - lcs);  // Can only be the prefix until the common suffix begins.
                --val.pos;
            }
            else
            {
                val.insBases = infix(alt, lcp, lcp + length(alt) - length(record.ref));
            }
            val.deltaType = DELTA_TYPE_INS;
        }
        else  // Simple deletion.
        {
            if (refBaseAfter)
            {
                val.insBases = prefix(record.ref, length(record.ref) - lcs);  // Can only be the prefix until the common suffix begins.
                --val.pos;
            }
            else
            {
                val.insBases = infix(record.ref, lcp, lcp + length(record.ref) - length(alt));
            }
            val.deltaType = DELTA_TYPE_DEL;
        }
    }
    return val;
}

inline String<BasicDeltaValue>
extractDeltaEvent(VcfRecord const & record, StringSet<CharString> const & altSet)
{
    String<BasicDeltaValue> deltas;

    // Extract delta event value for every ALT of the current record.
    for (auto & alt : altSet)
    {
        if (isSymbolicAllel(alt))  // Precise or imprecise structural variants
        {
            SEQAN_ASSERT_FAIL("Implement me!");
        }
        else if (isBreakendReplacement(alt))  // Complex breakpoints.
        {
            SEQAN_ASSERT_FAIL("Implement me!");
        }
        else  // Simple SNP or small InDel
        {
            appendValue(deltas, extractBasicDeltaValue(record, alt));
        }
    }
    return deltas;
}

inline String<BasicDeltaValue>
extractGenotypes(VcfRecord const & record)
{
    StringSet<CharString> altSet;
    strSplit(altSet, record.alt, EqualsChar<','>());

    // Extract genotypes per individual per haplotype.
    auto coverages = extractGenotypeInfos(record, altSet);
    auto deltaEvents = extractDeltaEvent(record, altSet);
    
    SEQAN_ASSERT_EQ(length(coverages), length(deltaEvents));
    
    auto itCov = begin(coverages, Standard());
    auto itDel = begin(deltaEvents, Standard());
    for (; itDel != end(deltaEvents, Standard()); ++itDel, ++itCov)
    {
        itDel->coverage = *itCov;
    }
    return deltaEvents;
}


// ----------------------------------------------------------------------------
// Function _readVcfRecords
// ----------------------------------------------------------------------------

template <typename TJst, typename TSize>
inline int
_readVcfRecords(VcfRecordTranslator<TJst> & delegate,
                VcfFileIn & vcfFile,
                TSize numSeq,
                ConverterOptions const & options)
{
//    typedef DeltaMap<TValue, TAlphabet1> TDeltaMap;
//    typedef typename DeltaValue<TDeltaMap, DeltaTypeIns>::Type TInsBuffer;
//    typedef typename DeltaValue<TDeltaMap, DeltaTypeDel>::Type TDel;
//    typedef VariantInfo<TDel, TDel, TInsBuffer> TVariantInfo;
//    typedef typename Value<TInsBuffer>::Type TAlphabet;


//#ifdef PRINT_DOTS
//    if (options.verbosity == 3)
    std::cerr << "\n# variants processed: ";
    static const unsigned DOT_SIZE = 10000;
    static const unsigned CHUNK_SIZE = 100000;
    unsigned pos = 0;
//#endif //PRINT_DOTS

////    unsigned numSeq = totalSequences;  // Total number of sequences.
//
//    String<bool, Packed<> > altCoverage;
//    resize(altCoverage, numSeq, false, Exact());   // Coverage of all sequences.

    if (options.includeReference)
        --numSeq;  // Remove reference from index because it cannot have an alternative site.
//    TVariantInfo _info;
//    _info.refPos = 0;
//    appendValue(_info.seqBuffer, TAlphabet());
//    _info.variantSize = 1;
//    _info.variantType = DELTA_TYPE_SNP;
//    appendValue(tmpInfo, _info);
//    appendValue(tmpCov, altCoverage);
    // Add dummy begin node.

    if (numSeq == 0)
        return 1;

    unsigned counter = 0;
    String<double> timeTable;
    resize(timeTable, 5, 0.0, Exact());

    CharString svFlag = "SVTYPE";
    VcfRecord vcfRecord;
    while(!atEnd(vcfFile))
    {
        clear(vcfRecord);
        ++counter;
        double timeRead = sysTime();
        readRecord(vcfRecord, vcfFile);

        timeTable[3] += sysTime() - timeRead;

//#ifdef PRINT_DOTS
        if (pos % CHUNK_SIZE == 0)
            std::cerr << " " << pos << " ";
//#endif //PRINT_DOTS

        // Only load records which are passed and have genotype information.
        toUpper(vcfRecord.filter);
        toUpper(vcfRecord.format);
        toUpper(vcfRecord.info);

        // Suppress SVs if enabled.
        if (options.suppressSVs)
        {

            if (std::search(begin(vcfRecord.format), end(vcfRecord.format), begin(svFlag), end(svFlag)) ==
                end(vcfRecord.format, Standard()))
                continue;
        }

//        if (vcfRecord.beginPos >= 38508)
//        {
//            std::cerr << "[LOG] POS " << vcfRecord.beginPos << " ALT " << vcfRecord.alt << " and ID " << vcfRecord.id << std::endl;
//        }

        // Check general filter options.
        if (vcfRecord.filter == "PASS" && startsWith(vcfRecord.format, "GT") &&
            !startsWith(vcfRecord.info, "IMPRECISE") && !(vcfRecord.qual != vcfRecord.qual))
        {
            String<BasicDeltaValue> basicDeltaValue = extractGenotypes(vcfRecord);

            if (empty(basicDeltaValue))  // No information extrated.
                continue;

            auto minHt = _min(length(options.haplotypes), length(basicDeltaValue));
            for (decltype(minHt) htNum = minHt; htNum < minHt; ++minHt)
            {
                double timeDelegate = sysTime();
                delegate(basicDeltaValue[htNum]);
                timeTable[2] += sysTime() - timeDelegate;
            }
//
//            // Extracts all alternative allels for the current loci.
//            String<TVariantInfo> variantInfoString;
//
//            double timeExtract = sysTime();
//            _extractVariants(variantInfoString, vcfRecord);
//            timeTable[0] += sysTime() - timeExtract;
//
//            // Coverages per alternative allel for the current loci.
//            String<String<bool, Packed<> > > altCoverageSet;
//            resize(altCoverageSet, length(variantInfoString), altCoverage, Exact());
//
//            // Read the combined genotype.
//            if (options.readGenotype)
//            {
//                for (unsigned seqId = 0; seqId < numSeq; ++seqId)
//                {
//                    // Process all haplotypes to generate a single genotype representing one value or not.
//                    int altPos = 0;
//                    for (unsigned htId = 0; htId < 2u; ++htId)  // NOTE(rmaerker): We only assume diploid.
//                    {
//                        double timeHaplotype = sysTime();
//                        int tmpAlt = _extractHaplotype(vcfRecord.genotypeInfos[seqId][htId << 1]);  // Extract the alt position for the current sequence.
//                        timeTable[1] += sysTime() - timeHaplotype;
//                        if(altPos == 0 && tmpAlt != 0)
//                        {
//                            altPos = tmpAlt;
//                            double timeAssignVec = sysTime();
//                            assignValue(altCoverageSet[tmpAlt - 1], seqId, true);
//                            timeTable[4] += sysTime() - timeAssignVec;
//                        }
//                        if (tmpAlt != 0 && tmpAlt != altPos)
//                            std::cerr << "WARNING: Different allels at the same loci in: " << vcfRecord.rID  << " at "  << vcfRecord.beginPos << " for " << seqId << "!" << std::endl;
//                    }
//                }
//            }
//            else
//            {
//                unsigned seqId = 0;
//                unsigned individualId = 0;
//                while (seqId < numSeq)
//                {
//                    // Iterate each haplotype.
//                    for (unsigned htId = 0; htId < length(options.haplotypes); ++htId, ++seqId)
//                    {
//                        // Extract the variant per individual per haplotype.
//                        double timeHaplotype = sysTime();
//                        int altPos = _extractHaplotype(vcfRecord.genotypeInfos[individualId][options.haplotypes[htId] * 2]);  // Extract the alt position for the current sequence.
//                        timeTable[1] += sysTime() - timeHaplotype;
//
//                        if (altPos != 0)
//                        {
//                            double timeAssignVec = sysTime();
//                            assignValue(altCoverageSet[--altPos], seqId, true);
//                            timeTable[4] += sysTime() - timeAssignVec;
//                        }
//                    }
//                    ++individualId;
//                }
//            }

//            std::cerr << "[LOG] Processed line: " << counter << std::endl;
//            if (counter == 1)
//            {
//                std::cerr << "[LOG] ________Print_START_________________________________________________" << std::endl;
//
//                for (unsigned j = 0; j < length(altCoverageSet); ++j)
//                {
//                    std::cerr << "[LOG] POS " << variantInfoString[j].refPos;
//                    switch(variantInfoString[j].variantType)
//                    {
//                        case DELTA_TYPE_SNP:
//                            std::cerr << " SNP " << variantInfoString[j].seqBuffer[0];
//                            break;
//                        case DELTA_TYPE_INS:
//                            std::cerr << " INS " << variantInfoString[j].seqBuffer;
//                            break;
//                        case DELTA_TYPE_DEL:
//                            std::cerr << " DEL " << variantInfoString[j].variantSize;
//                            break;
//                        default:
//                            SEQAN_THROW(-1);
//                    }
//                     std::cerr << " COV " << altCoverageSet[j] << std::endl;
//                }
//                std::cerr << "[LOG] ________Print_STOP__________________________________________________" << std::endl;
//            }
            // Iterates over the variants at this position.
        }
//#ifdef PRINT_DOTS
        if (++pos % DOT_SIZE == 0)
            std::cerr << ".";
//#endif //PRINT_DOTS
    }

    if (options.verbosity >= 2)
    {
        std::cerr << "\nNumber of records: " << counter << std::endl;
        std::cerr << "Time extract variant: " << timeTable[0] << " s." << std::endl;
        std::cerr << "Time extract haplotype: " << timeTable[1] << " s." << std::endl;
        std::cerr << "Time assign Vec: " << timeTable[4] << " s." << std::endl;
        std::cerr << "Time store data: " << timeTable[2] << " s." << std::endl;
        std::cerr << "Time read data: " << timeTable[3] << " s." << std::endl;
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function _convertVcfToGdf()
// ----------------------------------------------------------------------------

template <typename TReference, typename TRefAlphabet, typename TVarAlphabet>
int _convertVcfToGdf(GdfFileOut & gdfOut,
                     TReference & ref,
                     VcfFileIn & vcfFile,
                     ConverterOptions const & options,
                     TRefAlphabet const & /*refAlphabet*/,
                     TVarAlphabet const & /*snpAlphabet*/)
{
    typedef String<TRefAlphabet>        THost;
    typedef JournaledStringTree<THost>  TJst;

    TJst jst(ref, length(sampleNames(context(gdfOut))));

    VcfRecordTranslator<TJst> delegate(jst);

    double start = sysTime();
    _readVcfRecords(delegate, vcfFile, length(jst), options);
    if (options.verbosity > 1)
        std::cout << "Time for reading vcf: " << sysTime() - start << std::endl;

    _resolveConflicts(jst);

//    std::ofstream outputStream;
//    outputStream.open(toCString(converterOptions.outputFile), std::ios_base::out);
//    if (!outputStream.good())
//    {
//        std::cerr << "Cannot open file <"<< converterOptions.outputFile << "to write!";
//        return JSeqTools::FILE_READ_ERROR;
//    }

    start = sysTime();

    // Now we want to write the gdf file.

    // It would be nice to have such a file.


//    save(jst, gdfOut);  // Write save function. Export this later to develop.

    if (options.verbosity > 1)
    {
        std::cout << "Time for writing: " << sysTime() - start << std::endl;
        std::cout << "Number of converted nodes: " << size(jst);
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function adaptVcfToJournal()
// ----------------------------------------------------------------------------

template <typename TRefAlphabet, typename TVarAlphabet>
int adaptVcfToJournalData(ConverterOptions const & options,
                          TRefAlphabet const & /*refAlphabet*/,
                          TVarAlphabet const & /*varAlphabet*/)
{
    typedef String<TRefAlphabet, Alloc<> > TReference;

    if (options.verbosity > 0)
            std::cout << "Start converting vcf." << std::flush;

    // Open and read the sequence database.
    VcfFileIn vcfFile;

    if (!open(vcfFile, toCString(options.inputFile), OPEN_RDONLY))
    {
        std::stringstream msg;
        msg << "Could not open vcf file < " << options.inputFile << " >!";
        throw IOError(msg.str().c_str());
    }

    VcfHeader vcfHeader;
    readHeader(vcfHeader, vcfFile);

    unsigned numSamples = _min(length(sampleNames(context(vcfFile))), options.numIndividuals);
    // Resize to the maximum number of individuals that can be parsed.
    unsigned totalSequences = numSamples * length(options.haplotypes);

    GdfFileOut gdfOut;
    if (!open(gdfOut, toCString(options.outputFile), OPEN_WRONLY))
    {
        std::stringstream msg;
        msg << "Could not open gdf file < " << options.outputFile << " >!";
        throw IOError(msg.str().c_str());
    }
    // Load information for the gdfHeader.
    GdfHeader gdfHeader;

    // Prepare the sequence names for the output of the delta file.
    resize(sampleNames(context(gdfOut)), totalSequences, Exact());
    unsigned seqId = 0;
    for (unsigned individualId = 0; individualId < numSamples; ++individualId)
    {
        for (unsigned j = 0; j < length(options.haplotypes); ++j, ++seqId)
        {
            CharString fileName = getValue(sampleNames(context(vcfFile)), individualId);
            std::stringstream nameStream;
            nameStream << "_ht_" << options.haplotypes[j];
            append(fileName, nameStream.str());
            assignValue(sampleNames(context(gdfOut)), seqId, fileName);
        }
    }

   refresh(sampleNamesCache(context(gdfOut)));  // refresh the name store cache.

    // Load the reference infos.
    context(gdfOut).refUri = options.vcfReferenceFile;
    TReference reference;
    _loadContigs(context(gdfOut).refID, reference, toCString(context(gdfOut).refUri));

    // Add reference as file to the set of sequences.
    if (options.includeReference)
        appendValue(sampleNames(context(gdfOut)), context(gdfOut).refID);  // Last one.
    refresh(sampleNamesCache(context(gdfOut)));

    _convertVcfToGdf(gdfOut, reference, vcfFile, options, TRefAlphabet(), TVarAlphabet());

    if (options.verbosity > 0)
        std::cout << "\tDone!" << std::endl;

    return 0;
}

}

#endif // EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_ADAPT_VCF_H_