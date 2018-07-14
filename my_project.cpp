#include <seqan/bam_io.h>
#include "common.h"
#include "extension.h"
#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;
using namespace std;

myGlobalParameters params;


class myBamAlignmentRecord
{
public:
    CharString qName;               // QNAME
    uint16_t flag;                  // FLAG
    int32_t rID;                    // REF
    int32_t beginPos;               // POS
    uint8_t mapQ;                   // MAPQ mapping quality, 255 for */invalid
    uint16_t bin;                   // bin for indexing
    String<CigarElement<> > cigar;  // CIGAR string
    int32_t rNextId;                // RNEXT (0-based)
    int32_t pNext;                  // PNEXT (0-based)
    int32_t tLen;                   // TLEN
    CharString seq;                 // SEQ, as in SAM/BAM file.
    CharString qual;                // Quality string as in SAM (Phred).
    CharString tags;                // Tags, raw as in BAM.

    // Constants for marking pos, reference id and length members invalid (== 0/*).
    static int32_t const INVALID_POS = -1;
    static int32_t const INVALID_REFID = -1;
    static int32_t const INVALID_LEN = 0;
};

int main(int argc, char const ** argv)
{

    BamFileOut samFileOut(std::cout, Sam());
// BamFileOut samFileOut(context(bamFileIn), toCString(samFileOutName));
    
//     
//     typedef typename BamHeaderRecord::TTag    TTag;
/*    
@HD	VN:1.3	SO:coordinate
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
*/

    BamHeader header;
    BamHeaderRecord firstRecord;
    firstRecord.type = BAM_HEADER_FIRST;
    appendValue(firstRecord.tags, BamHeaderRecord::TTag("VN", "1.3"));
    appendValue(firstRecord.tags, BamHeaderRecord::TTag("SO", "coordinate"));
    appendValue(header, firstRecord);

    BamHeaderRecord seqRecord;
    seqRecord.type = BAM_HEADER_REFERENCE;
    appendValue(seqRecord.tags, BamHeaderRecord::TTag("SN", "ref"));
    appendValue(seqRecord.tags, BamHeaderRecord::TTag("LN", "45"));
    appendValue(seqRecord.tags, BamHeaderRecord::TTag("SN", "ref2"));
    appendValue(seqRecord.tags, BamHeaderRecord::TTag("LN", "40"));
    appendValue(header, seqRecord);
    
    writeHeader(samFileOut, header);
    
    vector<Pair<DnaString, Pair<uint32_t, uint32_t> > > occ;
    ;
    Pair<uint32_t, uint32_t> cord = Pair <uint32_t, uint32_t>(0, 32);
    occ.push_back(Pair < DnaString, Pair <uint32_t, uint32_t> > (DnaString("AAACATGCCCTCCCTTCCCCTGTGAGATTACTCGTGGGCAGGGCTTGTCCCTAGGGGGGCCCTGTCATTCAGTGGTTAACCCCGTGGTCCCAGGGTTTAC"), cord));
    cord = Pair <uint32_t, uint32_t>(0, 12);
    occ.push_back(Pair < DnaString, Pair <uint32_t, uint32_t> > (DnaString("TGTGCTAAGAATGTATTACTCTCTTTAGCAAATAACAGCCATAAATTTTGTTTCCAAAAATCTTTATCTTACTTGACAAGTGGGTGGAGACTGCTGAGCA"), cord));
    cord = Pair <uint32_t, uint32_t>(1, 62);
    occ.push_back(Pair < DnaString, Pair <uint32_t, uint32_t> > (DnaString("CCATAGTAGGTGCTCGGTAAATGTCAAAGATTTCCCAAATATTTGTATTGCTGCCAAGCTGAAGCATTTAGAAATGCATCAGTTTCATCTCCTCTTTGTG"), cord));
    cord = Pair <uint32_t, uint32_t>(1, 92);
    occ.push_back(Pair < DnaString, Pair <uint32_t, uint32_t> > (DnaString("AAACATGCCCTCCCTTCCCCTGTGAGATTACTCGTGGGCAGGGCTTGTCCCTAGGGGGGCCCTGTCATTCAGTGGTTAACCCCGTGGTCCCAGGGTTTAC"), cord));
    
    for(int i = 0; i < 4; ++i){
        
        BamAlignmentRecord record;
        record.qName = CharString("name" + to_string(i));
        record.rID = static_cast<uint32_t>(occ[i].i2.i1);
        record.beginPos = occ[i].i2.i2;
        record.seq = occ[i].i1;
        BamTagsDict tagsDict(record.tags);
        //getCigarString(...)
    
        setTagValue(tagsDict, "NM", 2);
        // => tags: "NM:i:2"
        setTagValue(tagsDict, "NH", 1);
        // => tags: "NM:i:2 NH:i:1"
        setTagValue(tagsDict, "NM", 3);
        // => tags: "NM:i:3 NH:i:1"
    
    writeRecord(samFileOut, record); 
    }
     
    
/*
    
    CharString bamFileInName = getAbsolutePath("demos/tutorial/sam_and_bam_io/example.sam");
    
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(bamFileInName)))
    {
        std::cerr << "ERROR: could not open input file " << bamFileInName << ".\n";
        return 1;
    }
    
//     BamFileOut samFileOut(context(bamFileIn), std::cout, Sam());

    BamHeader header;
    readHeader(header, bamFileIn);
    BamAlignmentRecord record;
    uint32_t count = 0;
    while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            cout << record.seq << endl;
            count += static_cast<int>(hasFlagUnmapped(record));
//             cout << std::string(record.cigar) << endl;
            
//             writeRecord(samFileOut, record); 
        }
    cout << "Count: " << count << endl;
  */
    
    
/*
    typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;
    
    TBamContext const & bamContext = context(bamFileIn);
    
    for(int i = 0; i< seqan::length(contigNames(bamContext)); ++i){
        cout << contigNames(bamContext)[i] << "\t" << contigLengths(bamContext)[i] << "\n";
    }*/
    
    
    
    /*
    try
    {
        BamHeader header;
        BamAlignmentRecord record;
        readHeader(header, bamFileIn);
        writeHeader(samFileOut, header);
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            writeRecord(samFileOut, record); 
        }
    }
    catch (ParseError const & e)
    {
        std::cerr << "ERROR: input header is badly formatted. " << e.what() << "\n";
    }
    catch (IOError const & e)
    {
        std::cerr << "ERROR: could not copy header. " << e.what() << "\n";
    }
    */
    
    
    /*
    CharString bamFileInName = getAbsolutePath("demos/tutorial/file_io_overview/example.bam");
    CharString samFileOutName = getAbsolutePath("demos/tutorial/file_io_overview/example.sam");

    // Open input BAM file.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(bamFileInName)))
    {
        std::cerr << "ERROR: could not open input file " << bamFileInName << ".\n";
        return 1;
    }

    // Open output SAM file.
//     BamFileOut samFileOut(context(bamFileIn), toCString(samFileOutName));
    BamFileOut samFileOut(context(bamFileIn), std::cout, Sam());

    // Copy header.
    BamHeader header;
    try
    {
        readHeader(header, bamFileIn);
        writeHeader(samFileOut, header);
    }
    catch (ParseError const & e)
    {
        std::cerr << "ERROR: input header is badly formatted. " << e.what() << "\n";
    }
    catch (IOError const & e)
    {
        std::cerr << "ERROR: could not copy header. " << e.what() << "\n";
    }

    // Copy all records.
    BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        try
        {
            readRecord(record, bamFileIn);
            writeRecord(samFileOut, record);
        }
        catch (ParseError const & e)
        {
            std::cerr << "ERROR: input record is badly formatted. " << e.what() << "\n";
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
    }*/

    return 0;
}

/*
Errors: 1   < AAACATGCCCTCCCTTCCCCTGTGAGATTACTCGTGGGCAGGGCTTGTCCCTAGGGGGGCCCTGTCATTCAGTGGTTAACCCCGTGGTCGCAGGGTTTAC , < 0 , 97437165 > >
AAACATGCCCTCCCTTCCCCTGTGAGATTACTCGTGGGCAGGGCTTGTCCCTAGGGGGGCCCTGTCATTCAGTGGTTAACCCCGTGGTCCCAGGGTTTAC
Errors: 1   < TGTGCTAAGAATGTATTAATCTCTTTAGCAAATAACAGCCATAAATTTTGTTTCCAAAAATCTTTATCTTACTTGACAAGTGGGTGGAGACTGCTGAGCA , < 0 , 97482311 > >
TGTGCTAAGAATGTATTACTCTCTTTAGCAAATAACAGCCATAAATTTTGTTTCCAAAAATCTTTATCTTACTTGACAAGTGGGTGGAGACTGCTGAGCA
Errors: 0   < CCATAGTAGGTGCTCGGTAAATGTCAAAGATTTCCCAAATATTTGTATTGCTGCCAAGCTGAAGCATTTAGAAATGCATCAGTTTCATCTCCTCTTTGTG , < 0 , 97547742 > >
CCATAGTAGGTGCTCGGTAAATGTCAAAGATTTCCCAAATATTTGTATTGCTGCCAAGCTGAAGCATTTAGAAATGCATCAGTTTCATCTCCTCTTTGTG
Errors: 2   < GTGCACATGCTCACACACATCCTCACACACATCCTTACACACCCTCACCCACATGCACTCACACACATGCACACACACTCCCTCACTCATGCACACATAC , < 0 , 97770616 > >
GTGCACATGCTCACACACATCCTCACACACTTCCTCACACACCCTCACCCACATGCACTCACACACATGCACACACACTCCCTCACTCATGCACACATAC*/

    /*
    myGlobalParameters myobject;
    params.flipdensity = 0.9;
    myobject.flipdensity = 0.8;
    cout << params.flipdensity << endl;
    cout << "Use Function: " << endl;
//     my::printAllP(params);
    params.print();
    cout << "Try function" << endl;
    cout << my::calcsomething(0.5, 0.5) << endl;
    */
 

