// vim:ts=8
#include <iostream>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <sstream>
#include <map>
#include <thread>
#include <chrono>
#include <gzstream.h>


#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>
#include "RGAssign.h"
#include "PutProgramInHeader.h"
#include "FastQParser.h"
//#include "SafeQueue.hpp"
#include "blockingconcurrentqueue.h"

#include "JSON.h"
#include "libgab.h"

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;
using namespace BamTools;

typedef struct { 

    ogzstream pairr1;
    ogzstream pairr2;
    ogzstream pairi1;
    ogzstream pairi2;

    ogzstream pairr1f;
    ogzstream pairr2f;
    ogzstream pairi1f;
    ogzstream pairi2f;


} fqwriters;

/****************************************/
/*                                      */
/*          STATIC VARIABLES            */
/*                                      */
/****************************************/


struct tallyForRG{    
    tallyForRG(void) : assigned(0), conflict(0), unknown(0), wrong(0) {}
    unsigned int assigned;
    unsigned int conflict;
    unsigned int unknown;
    unsigned int wrong;    
};

int qualOffset=33;

static string tagIndex1Seq  = "XI";
static string tagIndex1Qual = "YI";
static string tagIndex2Seq  = "XJ";
static string tagIndex2Qual = "YJ";

static double rgScoreCutoff  = 80 ;             // two HQ bases can mismatch
static double fracConflict   = 20 ;             // top shall be 100x more likely
static double wrongness      = 60 ;             // flag wrong if the wrong pair is a lot more likely
static int    mismatchesTrie = 2;
static int    maxErrorHits   = 20;
PrefixTree<string> * trieKnownString;
static string dashes = "--------------------------------";

map<string,string> p7id2seq;
map<string,string> p5id2seq;


// XXX really, a global?!
map<string,tallyForRG> namesMap; //a map name of RG to count of how many observed

// string getCWD(){
//    char temp[1000];
//    return ( getcwd(temp, 1000) ? string( temp ) : string("") );
// }

struct compareNameRG {
    bool operator() (pair<string,int> i,pair<string,int> j) {
        return (i.second>j.second);
    }
};

struct compareNameTally {
    bool operator() (pair<string,tallyForRG> i,pair<string,tallyForRG> j) { 
        return ( (i.second.assigned+i.second.unknown+i.second.conflict+i.second.wrong)
		 >
		 (j.second.assigned+j.second.unknown+j.second.conflict+j.second.wrong) );
    }
};


static string get_string_field( BamAlignment &al, const char* name ) 
{
    if(al.HasTag(name)) {
        char ttype;
        if( !al.GetTagType(name,ttype) ) {
            cerr << "Unable to get tag ("<<name<<") type for read  " << al.Name<<endl;
            exit(1);
        }
        if( ttype=='Z' || ttype=='H' ) {
            string tagInfo;
            if(!al.GetTag(name,tagInfo)){
                cerr << "Unable to edit " << name << " tag" << endl;
                exit(1);     
            }
            return tagInfo;
        } else if( ttype=='A' ) {
            int val;
            if( !al.GetTag(name,val) ){
                cerr << "Unable to edit " << name << " tag" << endl;
                exit(1);     
            }
            return string(1,val) ;
        }
    }
    return string() ;
}

static bool getIndices( const BamAlignment &al,string & index1,string & index1Q,string & index2,string & index2Q,bool suppressErrors=false){
    if(!al.GetTag(tagIndex1Seq,index1) ){ 	
        if (!suppressErrors) {
	    cerr << "Cannot retrieve "<<tagIndex1Seq<<" field  "<<al.Name << endl;
	    exit(1); 
        } else
            return false;
    }

    if(!al.GetTag(tagIndex1Qual,index1Q)){ 	
        if (!suppressErrors) {
	    cerr << "Cannot retrieve "<<tagIndex1Qual<<" field  "<<al.Name << endl;
	    exit(1); 
        } else
            return false;
    }

    if(!al.GetTag(tagIndex2Seq,index2)) {
        index2 ="";
        index2Q="";
    }
    else if(!al.GetTag(tagIndex2Qual,index2Q)) { //we got XJ, double indexed
        if (!suppressErrors) {
    	    cerr << "Cannot retrieve "<<tagIndex2Qual<<" field  "<<al.Name << endl;
    	    exit(1); 
        } else
            return false;
    }
    
    return true;
}


void updateRecord( BamAlignment &al, const rgAssignment &rg , bool failBAM)
{
    // get old ZQ field, remove "ICW"
    string zq = get_string_field(al, "ZQ");
    string::iterator p = zq.begin(), q = zq.begin(), e = zq.end() ;
    while( p != e ) {
        if( *p != 'I' && *p != 'C' && *p != 'W' ) 
        {
            *q = *p ;
            ++q ;
        }
        ++p ;
    }
    zq.erase( q, e ) ;

    string predictedGroup=rg.predictedGroup;
    //will overwrite the RG
    bool assigned=true;
    if( rg.predictedGroup.empty() ) {
	predictedGroup="unknown";
    }
    bool incrInTally=false;
    if( rg.conflict ){ zq += 'C' ; assigned=false; if(!incrInTally){ namesMap[predictedGroup].conflict++; incrInTally=true;}  }
    if( rg.wrong    ){ zq += 'W' ; assigned=false; if(!incrInTally){ namesMap[predictedGroup].wrong++;    incrInTally=true;}  }
    if( rg.unknown  ){ zq += 'I' ; assigned=false; if(!incrInTally){ namesMap[predictedGroup].unknown++;  incrInTally=true;}  }
    if(assigned){
	namesMap[predictedGroup].assigned++;
    }


    if( rg.predictedGroup.empty() ) {
        al.RemoveTag("RG");
        al.RemoveTag("Z0");
        al.RemoveTag("Z1");
        al.RemoveTag("Z2");
        zq += 'I';
	//namesMap[ "unknown" ] ++;
	//        al.EditTag("RG","Z","unknown");
    } else    {
	//namesMap[ predictedGroup ] ++;
        al.EditTag("RG","Z",rg.predictedGroup);

        al.EditTag("Z0","i",(int)round(-10 * rg.logLikelihoodScore));

        if( rg.logRatioTopToSecond <= 0 )
            al.EditTag("Z1","i",(int)round(-10 * rg.logRatioTopToSecond));
        else
            al.RemoveTag("Z1") ;



        if( (rg.topWrongToTopCorrect) >= 0)  //if the log odds ratio is greater than 0, the ratio is greater than 1 and the mispairing scenario is more likely than the correct pairing one
            al.EditTag("Z2","i",(int)round(10.0 * rg.topWrongToTopCorrect));
        else
            al.RemoveTag("Z2") ;

    }

    // store new ZQ field and set FailedQC flag if it isn't empty if failBAM
    if(failBAM){
	al.SetIsFailedQC( !zq.empty() ) ;
    
	if( zq.empty() ) 
	    al.RemoveTag("ZQ") ;
	else 
	    al.EditTag( "ZQ", "Z", zq ) ;
    }
}

inline bool containsNoNs(const string & sN){
    return (sN.find("N") == string::npos);
}

void initializeKnownIndices(PrefixTree<string> * trieKnownString,string configFile){
    //p7 300
    string line;
    ifstream myFile;
    string content="";
    myFile.open(configFile.c_str(), ios::in);

    if (myFile.is_open()){
	while ( getline (myFile,line)){
	    content+=line;
	}
	myFile.close();
    }else{
	cerr << "Unable to open config file "<<configFile<<endl;
	exit(1);
    }


    JSONValue *value = JSON::Parse(content.c_str());
    if (value == NULL){
	cerr<<"Failed to parse JSON file"<<endl;
	exit(1);
    }

    JSONObject root;
    root = value->AsObject();
    if(root.find(L"indices") == root.end()){
	cerr<<"Failed to parse JSON file, needs a indices field"<<endl;
	exit(1);
    }
    

    JSONValue *jsonIndices    = root.at(L"indices");
    JSONObject jsonIndicesObj = jsonIndices->AsObject();


    //p7 indices
    if(jsonIndicesObj.find(L"p7indices") == jsonIndicesObj.end()){
	cerr<<"Failed to parse JSON file, needs a p7indices field"<<endl;
	exit(1);
    }
        
    JSONValue *jsonIndicesp7    = jsonIndicesObj.at(L"p7indices");
    JSONObject jsonIndicesp7Obj = jsonIndicesp7->AsObject();
    JSONArray arrayp7 = jsonIndicesp7Obj[L"p7index"]->AsArray();
    for (unsigned int i = 0; i < arrayp7.size(); i++){
	JSONObject temp=	arrayp7[i]->AsObject();
	
	if(temp.find(L"seq") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a seq in the p7indices field"<<endl;
	    exit(1);
	}
	if(temp.find(L"id") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a id in the p7indices field"<<endl;
	    exit(1);
	}


	string tempSeq (temp[L"seq"]->AsString().begin(),
			temp[L"seq"]->AsString().end());
	string tempID  (temp[L"id"]->AsString().begin(),
			temp[L"id"]->AsString().end());
	trieKnownString->insertIntoTree( tempSeq.c_str() , "p7#"    +tempID);
	trieKnownString->insertIntoTree( tempSeq.c_str() , "p7REVC#"+tempID);
    }


    //p5 indices
    if(jsonIndicesObj.find(L"p5indices") == jsonIndicesObj.end()){
	cerr<<"Failed to parse JSON file, needs a p5indices field"<<endl;
	exit(1);
    }
        
    JSONValue *jsonIndicesp5    = jsonIndicesObj.at(L"p5indices");
    JSONObject jsonIndicesp5Obj = jsonIndicesp5->AsObject();
    JSONArray arrayp5 = jsonIndicesp5Obj[L"p5index"]->AsArray();
    for (unsigned int i = 0; i < arrayp5.size(); i++){
	JSONObject temp=	arrayp5[i]->AsObject();
	
	if(temp.find(L"seq") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a seq in the p5indices field"<<endl;
	    exit(1);
	}
	if(temp.find(L"id") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a id in the p5indices field"<<endl;
	    exit(1);
	}


	string tempSeq (temp[L"seq"]->AsString().begin(),
			temp[L"seq"]->AsString().end());
	string tempID  (temp[L"id"]->AsString().begin(),
			temp[L"id"]->AsString().end());

	trieKnownString->insertIntoTree( tempSeq.c_str() , "p5#"    +tempID);
	trieKnownString->insertIntoTree( tempSeq.c_str() , "p5REVC#"+tempID);
    }




    //Other sequences after the primming site
    string IS4="AGATCTC";
    trieKnownString->insertIntoTree( IS4.c_str() ,                    "IS4");
    trieKnownString->insertIntoTree( reverseComplement(IS4).c_str() , "REVC#IS4");


}


//todo
//add fastq

void readIndexNumbers(string configFile){


    //p7 300
    string line;
    ifstream myFile;
    string content="";
    myFile.open(configFile.c_str(), ios::in);

    if (myFile.is_open()){
	while ( getline (myFile,line)){
	    content+=line;
	}
	myFile.close();
    }else{
	cerr << "Unable to open config file "<<configFile<<endl;
	exit(1);
    }


    JSONValue *value = JSON::Parse(content.c_str());
    if (value == NULL){
	cerr<<"Failed to parse JSON file"<<endl;
	exit(1);
    }

    JSONObject root;
    root = value->AsObject();
    if(root.find(L"indices") == root.end()){
	cerr<<"Failed to parse JSON file, needs a indices field"<<endl;
	exit(1);
    }
    

    JSONValue *jsonIndices    = root.at(L"indices");
    JSONObject jsonIndicesObj = jsonIndices->AsObject();


    //p7 indices
    if(jsonIndicesObj.find(L"p7indices") == jsonIndicesObj.end()){
	cerr<<"Failed to parse JSON file, needs a p7indices field"<<endl;
	exit(1);
    }
        
    JSONValue *jsonIndicesp7    = jsonIndicesObj.at(L"p7indices");
    JSONObject jsonIndicesp7Obj = jsonIndicesp7->AsObject();
    JSONArray arrayp7 = jsonIndicesp7Obj[L"p7index"]->AsArray();
    for (unsigned int i = 0; i < arrayp7.size(); i++){
	JSONObject temp=	arrayp7[i]->AsObject();
	
	if(temp.find(L"seq") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a seq in the p7indices field"<<endl;
	    exit(1);
	}
	if(temp.find(L"id") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a id in the p7indices field"<<endl;
	    exit(1);
	}


	string tempSeq (temp[L"seq"]->AsString().begin(),
			temp[L"seq"]->AsString().end());
	string tempID  (temp[L"id"]->AsString().begin(),
			temp[L"id"]->AsString().end());
	//cerr<<"p7 #"<<tempID<<"#\t#"<<tempSeq<<"#"<<endl;
	if(p7id2seq.find(tempID) == p7id2seq.end()){
	    p7id2seq[tempID] = tempSeq;
	}else{
	    cerr<<"Error while parsing JSON file, ID: "<<tempID<<" was found twice"<<endl;
	    exit(1);
	}
	
    }


    //p5 indices
    if(jsonIndicesObj.find(L"p5indices") == jsonIndicesObj.end()){
	cerr<<"Failed to parse JSON file, needs a p5indices field"<<endl;
	exit(1);
    }
        
    JSONValue *jsonIndicesp5    = jsonIndicesObj.at(L"p5indices");
    JSONObject jsonIndicesp5Obj = jsonIndicesp5->AsObject();
    JSONArray arrayp5 = jsonIndicesp5Obj[L"p5index"]->AsArray();
    for (unsigned int i = 0; i < arrayp5.size(); i++){
	JSONObject temp=	arrayp5[i]->AsObject();
	
	if(temp.find(L"seq") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a seq in the p5indices field"<<endl;
	    exit(1);
	}
	if(temp.find(L"id") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a id in the p5indices field"<<endl;
	    exit(1);
	}


	string tempSeq (temp[L"seq"]->AsString().begin(),
			temp[L"seq"]->AsString().end());
	string tempID  (temp[L"id"]->AsString().begin(),
			temp[L"id"]->AsString().end());

	if(p5id2seq.find(tempID) == p5id2seq.end()){
	    p5id2seq[tempID] = tempSeq;
	}else{
	    cerr<<"Error while parsing JSON file, ID: "<<tempID<<" was found twice"<<endl;
	    exit(1);
	}


    }






}


void printUnfoundToFile(vector< pair<string,int> > * unfound,ofstream & fileError,const bool  printErrorJSONabsent){

    for(int i=0;i<min(int(unfound->size()),maxErrorHits);i++){	       
	//Searching in known strings
	vector<string> temp = allTokens((*unfound)[i].first,'#');
	vector<string> temp2;
	for(unsigned int j=0;j<temp.size();j++){
	    vector<string> * temp3=new vector<string>();
	    vector<string> * temp4=new vector<string>();
	    if(!printErrorJSONabsent){
		trieKnownString->searchMismatch(temp[j].c_str(),temp3,0);
		trieKnownString->searchMismatch( ( "N"+temp[j].substr(0, temp[j].size() -1)     ).c_str(),temp4,1);
		trieKnownString->searchMismatch( (     temp[j].substr(1, temp[j].size() -1)+"N" ).c_str(),temp4,1);

		//adding a tag before the shifted ones
		for(unsigned int k=0;k<temp4->size();k++)
		    (*temp4)[k]="SHFT#"+(*temp4)[k];
		temp3->insert( temp3->end(), temp4->begin(), temp4->end() );
	    }
	    if(temp3->size() == 0 && temp4->size() == 0){
		temp2.push_back( "?");
	    }else{
		temp2.push_back( vectorToString(*temp3,","));  
	    }


	    delete temp3;
	    delete temp4;

	    
	}
	
	fileError<< 
	    vectorToString( temp,"\t" )<<"\t"<<
	    (*unfound)[i].second<<"\t"<<
	    vectorToString( temp2,"\t" )<< endl;
    }
}

void check_thresholds( rgAssignment &rg ) {
    rg.unknown  =    (-10.0 * rg.logLikelihoodScore)   > rgScoreCutoff ;
    rg.conflict =  (( -10.0 * rg.logRatioTopToSecond)  < fracConflict) && (rg.logRatioTopToSecond < 0.0) ;
    rg.wrong    =  (   10.0 * rg.topWrongToTopCorrect) > wrongness ;
}

void processSingleEndReads( BamAlignment &al, BamWriter &writer, bool printError, map<string,int> &unknownSeq, map<string,int> &wrongSeq, map<string,int> &conflictSeq,bool failBAM)
{
    string index1;
    string index1Q;
    string index2;
    string index2Q;

    getIndices(al,index1,index1Q,index2,index2Q);

    rgAssignment rgReturn=assignReadGroup(index1,index1Q,index2,index2Q,rgScoreCutoff,fracConflict,mismatchesTrie,qualOffset);
    check_thresholds( rgReturn ) ;

    updateRecord(al,rgReturn,failBAM);
    writer.SaveAlignment(al);

    //record unresolved indices
    if(printError){
	string keyIndex;
	if(index2.empty()){
	    keyIndex=index1;
	}else{
	    keyIndex=index1+"#"+index2;
	}

	if( rgReturn.conflict ) conflictSeq[ keyIndex ] ++;
	if( rgReturn.unknown  ) unknownSeq [ keyIndex ] ++;
	if( rgReturn.wrong    ) wrongSeq   [ keyIndex ] ++;
    }

}

void processPairedEndReads( BamAlignment &al, BamAlignment &al2, BamWriter &writer, bool printError, map<string,int> &unknownSeq, map<string,int> &wrongSeq, map<string,int> &conflictSeq,bool failBAM)
{
    string index1;
    string index1Q;
    string index2;
    string index2Q;

    string sindex1;
    string sindex1Q;
    string sindex2;
    string sindex2Q;

    //retrieve indices
    getIndices(al,index1,index1Q,index2,index2Q);
    //check to see if the other indices are the same just for fun
    bool al2HasIndex=getIndices(al2,sindex1,sindex1Q,sindex2,sindex2Q,true);

    if(al2HasIndex && (index1 !=sindex1 )){cerr<<"Seq#1 has a different index 1 than seq #2, exiting "        <<al.Name<<" vs "<<al2.Name<< endl; exit(1);}
    if(al2HasIndex && (index1Q!=sindex1Q)){cerr<<"Seq#1 has a different index 1 quality than seq #2, exiting "<<al.Name<<" vs "<<al2.Name<< endl; exit(1);}
    if(al2HasIndex && (index2 !=sindex2 )){cerr<<"Seq#1 has a different index 2 than seq #2, exiting "        <<al.Name<<" vs "<<al2.Name<< endl; exit(1);}
    if(al2HasIndex && (index2Q!=sindex2Q)){cerr<<"Seq#1 has a different index 2 quality than seq #2, exiting "<<al.Name<<" vs "<<al2.Name<< endl; exit(1);}

    rgAssignment rgReturn = assignReadGroup(index1,index1Q,index2,index2Q,rgScoreCutoff,fracConflict,mismatchesTrie,qualOffset);
    check_thresholds( rgReturn ) ;

    updateRecord(al, rgReturn,failBAM);
    updateRecord(al2,rgReturn,failBAM);
    writer.SaveAlignment(al);
    writer.SaveAlignment(al2);

    // record unresolved indices
    if(printError){
        string keyIndex;
        if(index2.empty()){
            keyIndex=index1;
        }else{
            keyIndex=index1+"#"+index2;
        }

        if( rgReturn.conflict ) conflictSeq[ keyIndex ] += 1;
        if( rgReturn.unknown  ) unknownSeq [ keyIndex ] += 1;
        if( rgReturn.wrong    ) wrongSeq   [ keyIndex ] += 1;
    }
}


void checkFD(){
    int rofd  = returnOpenFileDescriptors();
    int rofdM = returnOpenFileDescriptorsMax();
    if(rofd == rofdM){
	cerr<<"WARNING: deML has detected that you have "<<rofd<<" open file descriptors out of a max. of "<<rofd<<endl;
	cerr<<"         therefore you have reached the maximum. If the information is correct, certain files might be empty"<<endl;
	cerr<<"         Either:"<<endl;
	cerr<<"           1) Use BAM as input/output"<<endl;
	cerr<<"           2) Check \"ulimit -n\" and put a higher number e.g. \"ulimit -n 1024\""<<endl;
	cerr<<"              If you already at the limit, increase the system limits"<<endl;
    }

}

struct FastQEntry{
    FastQEntry(kseq_t* kseq) : 
        name(kseq->name.s), 
        comment(kseq->comment.l ? kseq->comment.s : ""), 
        seq(kseq->seq.s),
        qual(kseq->qual.s)
    {}
    string name, comment, seq, qual;
};

class FastQFileReader{
public:
    FastQFileReader(const string& filename){
        fileHandle = gzopen(filename.c_str(), "r");
        kseq = kseq_init(fileHandle);
	tthread = std::thread(&FastQFileReader::run, this);
    }
    ~FastQFileReader(void){
        tthread.join();
	gzclose(fileHandle);
        delete kseq;
    }
    void run(void){
        FastQEntry* fqe;
        while(kseq_read(kseq) >=0){
	    fqe = new FastQEntry(kseq);
	    while(!queue.try_enqueue(fqe)){} //spin
        }
	queue.enqueue(0);
    }
    gzFile fileHandle;
    kseq_t* kseq;
    moodycamel::BlockingConcurrentQueue<FastQEntry*> queue;  
    std::thread tthread;
};

typedef std::tuple<FastQEntry*, FastQEntry*, FastQEntry*, FastQEntry*> FastQTuple;
typedef std::tuple<rgAssignment, string, string, string, string> ResultFastQTuple;

class Producer{
public:
    Producer(size_t max_queue_size, size_t n_clients, string r1f, string i1f, string r2f, string i2f) :
        queue(max_queue_size), n_clients(n_clients), forwardfq(r1f), index1fq(i1f), reversefq(r2f), index2fq(i2f)
    {}

    void operator()(void){
        FastQFileReader fqpf(forwardfq), fqpr(reversefq), fqpi1(index1fq), fqpi2(index2fq);
    	FastQEntry* ffo=0;
    	FastQEntry* rfo=0;
    	FastQEntry* i1fo=0;
    	FastQEntry* i2fo=0;

        fqpf.queue.wait_dequeue(ffo);
	while(ffo){
            fqpr.queue.wait_dequeue(rfo);
            fqpi1.queue.wait_dequeue(i1fo);
            fqpi2.queue.wait_dequeue(i2fo);
	    if(rfo==0||i1fo==0||i2fo==0){
	        cerr << "ERROR: Discrepancy between fastq files at record with first index " <<  ffo->name <<endl;
	        exit(1);
            }
            while(!queue.try_enqueue(FastQTuple(ffo, rfo, i1fo, i2fo))){} //stall while queue is full
	    fqpf.queue.wait_dequeue(ffo);
	}
        fqpr.queue.wait_dequeue(rfo);
        fqpi1.queue.wait_dequeue(i1fo);
        fqpi2.queue.wait_dequeue(i2fo);
	if(rfo!=0){
	    cerr << "ERROR: Discrepancy between fastq files at record with first index " <<  rfo->name <<endl;
	    exit(1);
        }
	if(i1fo!=0){
	    cerr << "ERROR: Discrepancy between fastq files at record with first index " <<  i1fo->name <<endl;
	    exit(1);
        }
	if(i2fo!=0){
	    cerr << "ERROR: Discrepancy between fastq files at record with first index " <<  i2fo->name <<endl;
	    exit(1);
        }
        for(size_t i=0; i<n_clients; ++i){
            //end markers
            queue.enqueue(FastQTuple(0,0,0,0));
        }
    }

    moodycamel::BlockingConcurrentQueue<FastQTuple> queue;  
//  SafeQueue<FastQTuple> queue;  
    size_t n_clients;
    string forwardfq, index1fq, reversefq, index2fq;
};

class FileThread{
public:
    FileThread(const string& prefixOut, const string& predictedGroup, bool hasRevBool, bool hasId2Bool, size_t n_producers) : queue((10 + 2 *  n_producers) * 32){
        //initialize writer
        rg2FqWriter = new fqwriters();

        string outpairr1   = prefixOut+"/"+predictedGroup+"_"+"r1.fq.gz";
        string outpairr2   = prefixOut+"/"+predictedGroup+"_"+"r2.fq.gz";
        string outpairi1   = prefixOut+"/"+predictedGroup+"_"+"i1.fq.gz";
        string outpairi2   = prefixOut+"/"+predictedGroup+"_"+"i2.fq.gz";

        string outpairr1f   = prefixOut+"/"+predictedGroup+"_"+"r1.fail.fq.gz";
        string outpairr2f   = prefixOut+"/"+predictedGroup+"_"+"r2.fail.fq.gz";
        string outpairi1f   = prefixOut+"/"+predictedGroup+"_"+"i1.fail.fq.gz";
        string outpairi2f   = prefixOut+"/"+predictedGroup+"_"+"i2.fail.fq.gz";

        rg2FqWriter->pairr1.open( outpairr1.c_str(),  ios::out);
        rg2FqWriter->pairr1f.open(outpairr1f.c_str(), ios::out);
	    	    
        rg2FqWriter->pairi1.open( outpairi1.c_str(),  ios::out);
        rg2FqWriter->pairi1f.open(outpairi1f.c_str(), ios::out);

        if(!rg2FqWriter->pairr1.good()){
            checkFD(); 
            cerr<<"Cannot write to file "<< outpairr1 <<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl;
            exit(1); 
        }
        if(!rg2FqWriter->pairr1f.good()){
            checkFD(); 
            cerr<<"Cannot write to file "<< outpairr1f<<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl; 
            exit(1); 
        }
	    	    
        if(!rg2FqWriter->pairi1.good()){
            checkFD(); 
            cerr<<"Cannot write to file "<< outpairi1 <<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl; 
            exit(1); 
        }
        if(!rg2FqWriter->pairi1f.good()){
            checkFD(); 
            cerr<<"Cannot write to file "<< outpairi1f<<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl;
            exit(1); 
        }

        if(hasId2Bool){
            rg2FqWriter->pairi2.open( outpairi2.c_str(),  ios::out);
            rg2FqWriter->pairi2f.open(outpairi2f.c_str(), ios::out);

            if(!rg2FqWriter->pairi2.good()){
                checkFD(); 
                cerr<<"Cannot write to file "<< outpairi2 <<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl;
                exit(1);
            }
            if(!rg2FqWriter->pairi2f.good()){
                checkFD(); cerr<<"Cannot write to file "<< outpairi2f<<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl; 
                exit(1); 
            }
        }

        if(hasRevBool){
            rg2FqWriter->pairr2.open( outpairr2.c_str(),  ios::out);
            rg2FqWriter->pairr2f.open(outpairr2f.c_str(), ios::out);

            if(!rg2FqWriter->pairr2.good()){
                checkFD(); 
                cerr<<"Cannot write to file "<<outpairr2 <<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl;
                exit(1);
            }
            if(!rg2FqWriter->pairr2f.good()){
                checkFD();
                cerr<<"Cannot write to file "<<outpairr2f<<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl; 
                exit(1); 
            }	   
        }
       
        //start thread
        tthread = std::thread(&FileThread::flush, this);
    }
    ~FileThread(void){
        ResultFastQTuple dummy;
        queue.enqueue(dummy);
        tthread.join(); //flush will end when it gets the dummy

        rg2FqWriter->pairr1.close( );
        rg2FqWriter->pairr1f.close();
	
        rg2FqWriter->pairi1.close( );
        rg2FqWriter->pairi1f.close();
	
        if(!rg2FqWriter->pairi2){
            rg2FqWriter->pairi2.close( );
        }
        if(!rg2FqWriter->pairi2f){
            rg2FqWriter->pairi2f.close();
        }
        if(!rg2FqWriter->pairr2){
            rg2FqWriter->pairr2.close( );
        }
        if(!rg2FqWriter->pairr2f){
            rg2FqWriter->pairr2f.close();	   
        }
        delete rg2FqWriter;
    }
    void enqueue(const ResultFastQTuple& rft){
        while(!queue.try_enqueue(rft)){
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }
    void flush(void){
        size_t count;
	bool done = false;
        while(!done){
            count = queue.wait_dequeue_bulk(buffer, buffer_size);
	    for(size_t i=0; i<count && !done; ++i){
                const ResultFastQTuple& rft(buffer[i]);
                if(get<1>(rft).empty()){
                    done = true;
                }
                else{
                    if(get<0>(rft).conflict  || get<0>(rft).wrong || get<0>(rft).unknown ){ //if the reads fail to meet thresholds
                        rg2FqWriter->pairr1f << get<1>(rft) <<endl;
                        rg2FqWriter->pairi1f << get<3>(rft) <<endl;
                        if(!get<4>(rft).empty()){
                            rg2FqWriter->pairi2f << get<4>(rft) <<endl;
                        }

                        if(!get<2>(rft).empty()){
                            rg2FqWriter->pairr2f << get<2>(rft) <<endl;
                        }
                    }else{
                        rg2FqWriter->pairr1 << get<1>(rft) <<endl;
                        rg2FqWriter->pairi1 << get<3>(rft) <<endl;
                        if(!get<4>(rft).empty()){
                            rg2FqWriter->pairi2 << get<4>(rft) <<endl;
                        }
	    
                        if(!get<2>(rft).empty()){
                            rg2FqWriter->pairr2 << get<2>(rft) <<endl;
                        }
                    }
                }
            }
        }
    }
    moodycamel::BlockingConcurrentQueue<ResultFastQTuple> queue;
    size_t buffer_size = 100;
    ResultFastQTuple buffer[100];
    fqwriters* rg2FqWriter;
    std::thread tthread;
};

class FileSorter{
public:
    FileSorter(string prefix, vector<string>& pgNames, bool hasRevBool, bool hasId2Bool, size_t num_worker_threads) {
        for(auto& pgName: pgNames){
            fileThreads[pgName] = new FileThread(prefix, pgName, hasRevBool, hasId2Bool, num_worker_threads);
        }
	fileThreads["unknown"] = new FileThread(prefix, "unknown", hasRevBool, hasId2Bool, num_worker_threads);
    }
    ~FileSorter(void){
        for(auto& kvp: fileThreads){
            delete kvp.second;
        }
    }
    void enqueue(const ResultFastQTuple& rft){
        const string& pg = get<0>(rft).predictedGroup;
        if(pg.empty()){
            fileThreads.at("unknown")->enqueue(rft);
        }
	else{
            fileThreads.at(pg)->enqueue(rft);
        }
    }
private:
    map<string,FileThread*> fileThreads;
};

string formatFastQ(FastQEntry* kseq){
    string ret = "";
    ret += kseq->name;
    if(!kseq->comment.empty()){
    	ret += " ";
    	ret += kseq->comment;
    }
    ret += "\n";
    ret += kseq->seq;
    ret += "\n+\n";
    ret += kseq->qual;
    ret += "\n";
    return ret;
}

class FastQWorker{
public:
    FastQWorker(moodycamel::BlockingConcurrentQueue<FastQTuple>& in_queue, FileSorter& out_queue, bool printSummary, bool printError) :
        in_queue(in_queue), out_queue(out_queue), 
        printSummary(printSummary), printError(printError)
    {}
    
    void operator()(void){
        FastQEntry* ffo=0;
        FastQEntry* rfo=0;
        FastQEntry* i1fo=0;
        FastQEntry* i2fo=0;
	string ffo_s, rfo_s, i1fo_s, i2fo_s;
	string index1s, index1q, index2s, index2q;
	FastQTuple ft;
	
        while(true){
            in_queue.wait_dequeue(ft);
            std::tie(ffo, rfo, i1fo, i2fo) = ft;
            if(ffo==0 && rfo==0 && i1fo==0 && i2fo==0)
                break;

            string deffs = ffo->name;
            if(strEndsWith(deffs,  "/1")){
                deffs=deffs.substr(0,deffs.size()-2);
            }

            string defi1s = i1fo->name;

            if( (deffs != defi1s ) ){
               cerr << "ERROR: Discrepancy between fastq files, different names with first index " <<deffs <<" and "<< defi1s <<endl;
               exit(1);
            }

            index1s = i1fo->seq;
            index1q = i1fo->qual;

            if(i2fo){
                string defi2s = i2fo->name;
                if( (deffs != defi2s ) ){
                    cerr << "ERROR: Discrepancy between fastq files, different names with second index " <<deffs <<" and "<< defi2s <<endl;
                    exit(1);
                }

                index2s =  i2fo->seq;
                index2q =  i2fo->qual;
            }

	    if(rfo){
                string defrs = rfo->name;
	    
                if(strEndsWith(defrs,  "/2")){
                    defrs=defrs.substr(0,defrs.size()-2);
                }

                if( (deffs != defrs) ){
                    cerr << "ERROR: Discrepency between fastq files, different names " <<deffs <<" and "<< defrs <<endl;
                    exit(1);
                }
            }

            rgAssignment rgReturn = assignReadGroup(index1s,index1q,index2s,index2q,rgScoreCutoff,fracConflict,mismatchesTrie,qualOffset);
            check_thresholds( rgReturn ) ;
            string predictedGroup;

            //////////////////////////
            //BEGIN update counters //
            //////////////////////////
            if(printSummary){
                string predictedGroup=rgReturn.predictedGroup;
                //will overwrite the RG
                bool assigned=true;
                if( rgReturn.predictedGroup.empty() ) {
                    predictedGroup="unknown";
                }
                bool incrInTally=false;
                if( rgReturn.conflict ){
                    assigned=false; 
                    if(!incrInTally){
                        namesMap[predictedGroup].conflict++; 
                        incrInTally=true;
                    }
                }
                if( rgReturn.wrong    ){
                    assigned=false; 
                    if(!incrInTally){ 
                        namesMap[predictedGroup].wrong++;
                        incrInTally=true;
                    }  
                }
                if( rgReturn.unknown  ){
                    assigned=false; 
                    if(!incrInTally){
                        namesMap[predictedGroup].unknown++;
                        incrInTally=true;
                    }
                }
                if(assigned){
                    namesMap[predictedGroup].assigned++;
                }
            }
	
            if(printError){
                string keyIndex;
                if(!i2fo){
                    keyIndex=index1s;
                }else{
                    keyIndex=index1s+"#"+index2s;
                }
                if( rgReturn.conflict ) conflictSeq[ keyIndex ] += 1;
                if( rgReturn.unknown  ) unknownSeq [ keyIndex ] += 1;
                if( rgReturn.wrong    ) wrongSeq   [ keyIndex ] += 1;
            }
            //////////////////////////
            //  END update counters //
            //////////////////////////
	    ffo_s = formatFastQ(ffo);
	    if(rfo){rfo_s = formatFastQ(rfo);} else {rfo_s="";}
	    i1fo_s = formatFastQ(i1fo);
	    if(i2fo){i2fo_s = formatFastQ(i2fo);} else {i2fo_s="";}

      	    out_queue.enqueue(ResultFastQTuple(rgReturn, ffo_s, rfo_s, i1fo_s, i2fo_s));
	    delete ffo;
	    delete i1fo;
	    if(rfo) delete rfo;
	    if(i2fo) delete i2fo;
	}
    }

    moodycamel::BlockingConcurrentQueue<FastQTuple>& in_queue;  
//    SafeQueue<FastQTuple>& in_queue;
//    SafeQueue<ResultFastQTuple>& out_queue;
    FileSorter& out_queue;  
    map<string,tallyForRG> namesMap; //a map name of RG to count of how many observed
    bool printSummary, printError;
    map<string,int> conflictSeq, unknownSeq, wrongSeq;
};

void processFastq(string           forwardfq,
		  string           reversefq,
		  string           index1fq,
		  string           index2fq,
		  string           prefixOut,
		  map<string,int> &unknownSeq, 
		  map<string,int> &wrongSeq, 
		  map<string,int> &conflictSeq,
		  const bool       printSummary,
		  const bool       printError,
		  size_t           num_worker_threads,
		  size_t           queue_size,
                  vector<string>&   pgNames){

    size_t max_queue_size = num_worker_threads * queue_size;

    //initialize thread objects
    Producer producer(max_queue_size, num_worker_threads, forwardfq, index1fq, reversefq, index2fq);
    FileSorter fileSorter(prefixOut, pgNames, !reversefq.empty(), !index2fq.empty(), num_worker_threads);
    vector<FastQWorker> workers(
        num_worker_threads,
        FastQWorker(producer.queue, fileSorter, printSummary, printError)
    );

    //start threads
    std::thread producer_thread(std::ref(producer));
    vector<std::thread> worker_threads;
    for(auto& worker: workers){
        worker_threads.push_back(std::thread(std::ref(worker)));
    }

    //join threads
    producer_thread.join();
    for(auto& worker_thread: worker_threads){
        worker_thread.join();
    }
  
    //collect stats
    for(auto& worker: workers){
        for(auto const&x: worker.unknownSeq){
            unknownSeq[x.first] += x.second;
        }
        for(auto const&x: worker.wrongSeq){
            wrongSeq[x.first] += x.second;
        }
        for(auto const&x: worker.conflictSeq){
            conflictSeq[x.first] += x.second;
        }
        for(auto const&x: worker.namesMap){
            namesMap[x.first].assigned += x.second.assigned;
            namesMap[x.first].conflict += x.second.conflict;
            namesMap[x.first].unknown += x.second.unknown;
            namesMap[x.first].wrong += x.second.wrong;
        }
    }
    //fileSorter flushes here
}



int main (int argc, char *argv[]) {

    BamReader reader;
    BamWriter writer;

    string bamFile;
    string bamFileOUT="";

    string index="";
    string outfile;
    bool   printSummary=false;
    string filenameSummary;

    bool   printError=false;
    bool   printErrorJSONabsent=false;

    string filenameError;

    ofstream ratioValuesOS;
    ofstream rgqualOS;

    bool ratioValuesFlag = false;
    bool rgqualFlag      = false;
    bool shiftByOne      = false;
    bool failBAM         = true;

    bool useFastq=false;
    string forwardfq;
    string reversefq;
    string index1fq;
    string index2fq;

    size_t num_worker_threads = 1;
    size_t queue_size = 10;

    bool produceUnCompressedBAM=false; 

    const string usage=string(string(argv[0])+
			      " [options] BAMfile"+"\n\n"+

			      "\t"+"Mandatory parameters:"+"\n"+
			      "\t\t"+"-i"+","+"--index"+"\t[index]"+"\t\t\t"+"File describing index sequences used"+"\n"+
			      "\t\t"+"-o"+","+"--outfile"+"\t[outfile]"+"\t\t"+"Specify output file"+"\n\n"+

			      "\t"+"BAM input options:"+"\n"+
			      "\t\t"+"--bamtags"+"\t[Idx1Seq],[Idx1Qual],[Idx2Seq],[Idx2Qual]"+"\t"+"BAM Tags containing barcodes\n"+

			      "\t"+"Fastq input/output (optional):"+"\n"+
			      "\t\t"+"You can specify fastq as input/output, in which case the -o option will be"+"\n"+
			      "\t\t"+"treated as an output prefix"+"\n"+
			      
			      "\t\t"+"-f"  +"\t[forward fastq]"+"\t\t\t"+""+"Forward reads in fastq\n"+
			      "\t\t"+"-r"  +"\t[reverse fastq]"+"\t\t\t"+""+"Reverse reads in fastq\n"+
			      "\t\t"+"-if1"  +"\t[index1 file fastq]"+"\t\t"+""+"First index sequences in fastq\n"+
			      "\t\t"+"-if2"  +"\t[index2 file fastq]"+"\t\t"+""+"Second index sequences in fastq\n"+
			      
			      "\tCutoffs options:"+"\n"
			      "\t\t"+"--phred64"  +"\t"+"\t\t\t"+""+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+

			      "\t\t"+"--rgqual"  +"\t[quality]"+"\t\t"+""+"Worst quality before flagging as unknown ["+stringify(rgScoreCutoff)+"]\n"+
			      "\t\t"+"--fracconf"+"\t[quality]"+"\t\t"+""+"Maximum quality difference before being considered a conflict ["+stringify(fracConflict)+"] \n"+
			      "\t\t"+"--wrongness"+"\t[quality]"+"\t\t"+""+"Mininum quality difference to flag as wrongly paired ["+stringify(wrongness)+"] \n"+
			      "\t\t"+"--mm"+"\t\t[mismatches]"+"\t\t"+""+"Maximum # of tolerated mismatches ["+stringify(mismatchesTrie)+"] \n"+

			      "\n\tRG assignment options:"+"\n"+
			      "\t\t"+"" +""+"--shift"+"\t"+"\t\t\t\t"+"Try shifting the index right by one at the cost of a mismatch"+"\n"+

			      
			      //"\t"+"\tOptional:"+"\n"+
			      "\n\tOutput options:"+"\n"+
			      "\t\t"+"--maxerr"+"\t[max err]"+"\t\t"+""+"Print  # wrongly of assigned RG in the error log (--error) ["+stringify(maxErrorHits)+"] \n"+
                              "\t\t"+"-u" +"\t\t\t\t\t"           +"Produce uncompressed bam (good for pipe)"+"\n"+ 
			      "\t\t"+"-s"+","+"--summary"+"\t[summary file]"+"\t\t"+"Summarize the RG tally in this file"+"\n"+
			      "\t\t"+"-e"+","+"--error"  +"\t[error file]"+"\t\t"+"Summarize the indices that were not assigned to a RG"+"\n"+
			      "\t\t"+""+""+"--rgval"  +"\t[file]"+"\t\t\t\t"+"Write the rg qualities as a binary file"+"\n"+
			      "\t\t"+""+""+"--ratio"   +"\t\t[file]"+"\t\t\t"+"Write the likelihood ratios as a binary file"+"\n"+
			      "\t\t"+""+""+"--nofail"   +"\t\t\t"+"\t"+"Do not set the QC fail in the output BAM file"+"\n"+
			      "\n\tThreadings options:"+"\n"+
                              "\t\t"+"-t" +"\t\t\t\t\t"           +"Number of threads"+"\n"+ 
                              "\t\t"+"-q" +"\t\t\t\t\t"           +"Queue size"+"\n"

			      );
			      

    if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
    	cout<<"Usage:"<<endl;
    	cout<<""<<endl;
    	cout<<usage<<endl;
    	return 1;
    }

    int lastIndexArgc = argc-1;
    for(int i=1;i<(lastIndexArgc);i++){ //all but the last
	
	if(string(argv[i]) == "--nofail"  ){
            failBAM=false;
            continue;
        }

	if(strcmp(argv[i],"-f") == 0 ){
	    forwardfq = string(argv[i+1]);
	    useFastq=true;
	    lastIndexArgc = argc;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-r") == 0 ){
	    reversefq = string(argv[i+1]);
	    useFastq=true;
	    lastIndexArgc = argc;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--bamtags") == 0 ){
	    istringstream iss(argv[i+1]);
	    string tok;
	    int p=0;
	    while(getline(iss, tok, ',')) {
	        switch (p) {
	            case 0: tagIndex1Seq = tok; break;
	            case 1: tagIndex1Qual = tok; break;
	            case 2: tagIndex2Seq = tok; break;
	            case 3: tagIndex2Qual = tok; break;
	            default:
	                cerr<<"too many commas in "<<argv[i+1]<<endl;
	                exit(1);
	        }
	        ++p;
	    }
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-if1") == 0 ){
	    index1fq = string(argv[i+1]);
	    useFastq=true;
	    lastIndexArgc = argc;
	    i++;
	    continue;
	}

	if(string(argv[i]) == "--phred64"  ){
            qualOffset=64;
            continue;
        }

	if(strcmp(argv[i],"-if2") == 0 ){
	    index2fq = string(argv[i+1]);
	    useFastq=true;
	    lastIndexArgc = argc;
	    i++;
	    continue;
	}
	

	if(strcmp(argv[i],"--shift") == 0 ){
	    shiftByOne      = true;
	    continue;
	}

	if(strcmp(argv[i],"--rgval") == 0 ){
	    string temp =string(argv[i+1]);
	    rgqualOS.open(temp.c_str(), ios::out | ios::binary);
	    rgqualFlag      = true;
	    if (!rgqualOS){
		cerr<<"Cannot print to file "<<temp<<endl;
		exit(1);
	    }
	    setFileForRGQual(&rgqualOS);
	    i++;
	    continue;
	}


	if(strcmp(argv[i],"-u") == 0  ){ 
	    produceUnCompressedBAM=true; 
	    continue; 
	} 
	

	if(strcmp(argv[i],"--ratio") == 0 ){
	    string temp =string(argv[i+1]);
	    ratioValuesOS.open(temp.c_str(), ios::out | ios::binary);
	    ratioValuesFlag = true;

	    if (!ratioValuesOS){
		cerr<<"Cannot print to file "<<temp<<endl;
		exit(1);
	    }
	    setFileForRatio(&ratioValuesOS);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-e") == 0 || strcmp(argv[i],"--error") == 0 ){
	    printError=true;
	    filenameError =string(argv[i+1]);
	    i++;
	    continue;
	}



	if(strcmp(argv[i],"-s") == 0 || strcmp(argv[i],"--summary") == 0 ){
	    printSummary=true;
	    filenameSummary =string(argv[i+1]);
	    i++;
	    continue;
	}





	if(strcmp(argv[i],"--maxerr") == 0 ){
	    maxErrorHits =destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}


	if(strcmp(argv[i],"--rgqual") == 0 ){
	    rgScoreCutoff =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}



	if(strcmp(argv[i],"--fracconf") == 0 ){
	    fracConflict =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--wrongness") == 0 ){
	    wrongness =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--mm") == 0 ){
	    mismatchesTrie =destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}



	if(strcmp(argv[i],"-i") == 0 || strcmp(argv[i],"--index") == 0 ){
	    index =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--outfile") == 0 ){
	    outfile=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-t") == 0){
            num_worker_threads = destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-q") == 0){
            queue_size = destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;             
    }



    if(useFastq){

	if(forwardfq.empty() ){
	    cerr<<"If fastq is used, the forward read must be specified"<<endl;
	    return 1;             
	}

	if(index1fq.empty() ){
	    cerr<<"If fastq is used, the forward read must be specified"<<endl;
	    return 1;             
	}

	if(produceUnCompressedBAM){
	    cerr<<"If fastq is used, cannot use -u flag"<<endl;
	    return 1;             
	}
	
    }
    
    
    
    if(outfile.size() == 0){
	cerr<<"The field -o is mandatory exiting"<<endl;
	return 1;             
    }

    if(index.size() == 0){
	cerr<<"The field -i is mandatory exiting"<<endl;
	return 1;             
    }


    if(printError){
	trieKnownString = new PrefixTree<string>();

	string configFile = "/../../webForm/config.json";

	if(!isFile( getCWD(argv[0])+configFile )){
	    configFile = "/../webForm/config.json";
	    if(!isFile( getCWD(argv[0])+configFile )){
		cerr<<"ERROR: file "<<(getCWD(argv[0])+"/../../webForm/config.json")<<" or "<<(getCWD(argv[0])+"/../webForm/config.json")<<" were not found, this file is required with the -e or --error option"<<endl;
		//return 1;             
	    }
	    printErrorJSONabsent=true;
	}else{
	    initializeKnownIndices(trieKnownString,getCWD(argv[0])+configFile);
	}
	// //debug
	// vector<string> * temp3=new vector<string>();
	// vector<string> * temp4=new vector<string>();
	// string t="AAGGTCT";
	// trieKnownString->searchMismatch(t.c_str(),temp3,0);
	// trieKnownString->searchMismatch( ( "N"+t.substr(0, t.size() -1)     ).c_str(),temp4,1);
	// trieKnownString->searchMismatch( (     t.substr(1, t.size() -1)+"N" ).c_str(),temp4,1);
	// cout<<"t3 "<<vectorToString(*temp3)<<endl;
	// cout<<"t4 "<<vectorToString(*temp4)<<endl;

	// return 1;
	// //end debug
    }

    bamFile=argv[argc-1];
    ifstream myIndexFile;


    //
    //  BEGIN : reading index file
    //    
    string indexStringFile="";
    bool firstLine=true;
    bool fileContainsSeq=true;
    bool fileContainsSeqDouble=false;

    myIndexFile.open(index.c_str(), ios::in);
    if (myIndexFile.is_open()){
	string line;
	while ( getline (myIndexFile,line)){
	    indexStringFile+=line+"\n";

	    if(firstLine){
		firstLine=false;
	    }else{
		vector<string> tempf = splitWS(line);

		if(tempf.size() == 2){ //contains 2 fields, single index
		    fileContainsSeqDouble=false;

		    if(isStringNatNumber(    tempf[0]) || 
		       isStringNatNumber("t"+tempf[0])){ //for truseq
			fileContainsSeq=false;
		    }

		}else{
		    if(tempf.size() == 3){ //contains 3 fields, double index
			fileContainsSeqDouble=true;

			if(isStringNatNumber(    tempf[0]) || 
			   isStringNatNumber("t"+tempf[0])){ //for truseq
			    fileContainsSeq=false;

			    if( isStringNatNumber(    tempf[1]) || 
			        isStringNatNumber("t"+tempf[1]) ){ //for truseq
				//fine
			    }else{
				cerr << "Error: line "<<line<<" cannot have a mixture of index sequence and index numbers"<<endl;
				exit(1);				
			    }

			}

		    }else{
			cerr << "Error: line "<<line<<" does not have 2 or 3 fields"<<endl;
			exit(1);
		    }

		}

	    }
	}
	myIndexFile.close();
    }else{ 
	cerr << "Unable to open file "<<index<<endl;
	exit(1);
    }
    //
    //  END : reading index file
    //    


    //
    //  BEGIN : Processing index with numbers
    //
    if(!fileContainsSeq){


	string configFile = "/../../webForm/config.json";

	if(!isFile( getCWD(argv[0])+configFile )){
	    configFile = "/../webForm/config.json";
	    if(!isFile( getCWD(argv[0])+configFile )){
		cerr<<"ERROR: file "<<(getCWD(argv[0])+"/../../webForm/config.json")<<" or "<<(getCWD(argv[0])+"/../webForm/config.json")<<" were not found, this file is required with the -e or --error option"<<endl;
		return 1;             
	    }
	}
	//initializeKnownIndices(trieKnownString,getCWD(argv[0])+configFile);
	
	readIndexNumbers( getCWD(argv[0])+configFile);	

	firstLine=true;
	vector<string> allLinesIndexFile = allTokens(indexStringFile,'\n');
	string tempIndexStringFile="";
	//while ( getline (myFile,line)){
	for(unsigned int i=0;i<allLinesIndexFile.size();i++){
	    string line = allLinesIndexFile[i];
	    if(line.empty())
		continue;
	    vector<string> tempf = splitWS(line);

	    if( ((tempf.size() == 2) && !fileContainsSeqDouble) || 
		((tempf.size() == 3) && fileContainsSeqDouble) ){
	    }else{
		cerr << "Error: line "<<line<<" does not have 2 or 3 fields"<<endl;
		exit(1);
	    }

	    if(firstLine){
		tempIndexStringFile += line;
		firstLine=false;
	    }else{

		if(p7id2seq.find(tempf[0]) != p7id2seq.end()){
		    tempIndexStringFile +=  p7id2seq[ tempf[0] ] ;
		}else{
		    cerr<<"Error p7 ID: #"<<tempf[0]<<"# not found"<<endl;
		    exit(1);
		}
	    
		if(fileContainsSeqDouble){

		    if(p5id2seq.find(tempf[1]) != p5id2seq.end()){
			tempIndexStringFile +=  "\t"+p5id2seq[ tempf[1] ] ;
		    }else{
			cerr<<"Error p5 ID: #"<<tempf[1]<<"# not found"<<endl;
			//cerr<<"Error while parsing JSON file, ID: "<<tempID<<" was found twice"<<endl;
			exit(1);
		    }
		
		    tempIndexStringFile +=  "\t"+tempf[2]  ;
		}else{
		    tempIndexStringFile +=  "\t"+tempf[1]  ;
		}
	    }
	    
	    tempIndexStringFile +=  "\n";
	    //cerr<<"temp "<<tempIndexStringFile<<endl;
	}//for all lines
	
	indexStringFile=tempIndexStringFile;
	cerr<<"Using the following as index sequences"<<endl;
	cerr<<"----------------------------------"<<endl;
	cerr<<indexStringFile<<endl;
	cerr<<"----------------------------------"<<endl;
    }//end file contains seq
    //
    //  END : Processing index with numbers
    //



    map<string,string> rgs = readIndexFile(indexStringFile,mismatchesTrie,shiftByOne);
    map<string,int> unknownSeq;
    map<string,int> wrongSeq;
    map<string,int> conflictSeq;
    
    if(useFastq){
        vector<string> rgNames;
        for(auto& kvp: rgs){
            rgNames.push_back(kvp.first);
        }

	processFastq(forwardfq,
		     reversefq,
		     index1fq,
		     index2fq,
		     outfile,
		     unknownSeq,
		     wrongSeq,
		     conflictSeq,
		     printSummary,
		     printError,
		     num_worker_threads,
		     queue_size,
                     rgNames);

    }else{


	if ( !reader.Open(bamFile) ) {
	    cerr << "Could not open input BAM file  "<<bamFile << endl;
	    return 1;
	}

	SamHeader  myHeader=reader.GetHeader();
	SamProgram sp;
   
	string pID          = "deML";   
	string pName        = "deML";   
	string pCommandLine = "";
	for(int i=0;i<(argc);i++){
	    pCommandLine += (string(argv[i])+" ");
	}

	putProgramInHeader(&myHeader,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),".."));

	SamReadGroupDictionary  srgd;
	map<string,string>::const_iterator itRG;   
	for ( itRG=rgs.begin(); itRG != rgs.end(); itRG++ ){
	    SamReadGroup srg ( itRG->first );	
	    srg.Description  = itRG->second; //description read in index file
	    srgd.Add( srg );       	

	    namesMap[ itRG->first ].assigned =0;
	    namesMap[ itRG->first ].unknown  =0;
	    namesMap[ itRG->first ].conflict =0;
	    namesMap[ itRG->first ].wrong    =0;


	    if(itRG->first == "conflict" || itRG->first == "unknown" || itRG->first == "wrong" ){
		cerr<<"ERROR: The RG names cannot contain the words: \"conflict\" or \"unknown\" or \"wrong\""<<endl;
		return 1;
	    }
	}

	namesMap[ "unknown" ].assigned=0;
	namesMap[ "unknown" ].unknown=0;
	namesMap[ "unknown" ].conflict=0;
	namesMap[ "unknown" ].wrong=0;



	myHeader.ReadGroups=srgd;
	if(produceUnCompressedBAM)  
	    writer.SetCompressionMode(BamWriter::Uncompressed); 

	if( !writer.Open(outfile,myHeader,reader.GetReferenceData() ) ) {
	    cerr << "Could not open output BAM file  "<<outfile << endl;
	    return 1;	
	}


	BamAlignment al;
	BamAlignment al2;

	// damn, this logic is convoluted...
	while( reader.GetNextAlignment(al) ) {
	    while(1) {
		if( !reader.GetNextAlignment(al2) ) {
		    // EOF, process the one leftover record
		    processSingleEndReads(al,writer,printError,unknownSeq,wrongSeq,conflictSeq,failBAM);
		    break; 
		}
		// If it's paired, both should have the same index, and we
		// save some work.  Since the reads are probably not
		// ordered, check the names first
		if( al.IsPaired() && al.Name == al2.Name ) {
		    processPairedEndReads(al,al2,writer,printError,unknownSeq,wrongSeq,conflictSeq,failBAM);
		    break ;
		} else {
		    // no match, treat one(!) separately
		    processSingleEndReads(al ,writer,printError,unknownSeq,wrongSeq,conflictSeq,failBAM);
		    swap(al,al2) ;
		}
	    }
	}

	reader.Close();
	writer.Close();






    }//bam
    //cleaning up

    //Print summary of RG assignment
    if(printSummary){
	map<string,tallyForRG>::iterator it;   
	unsigned int totalRG=0;	
	unsigned int totalAssignRG=0;	

	vector< pair<string,tallyForRG> > toprintVec;
	for ( it=namesMap.begin() ; it != namesMap.end(); it++ ){
	    toprintVec.push_back(  make_pair( it->first , it->second ) );
	    totalRG+=it->second.assigned+it->second.unknown+it->second.conflict+it->second.wrong;
	}

	sort (toprintVec.begin(),   toprintVec.end(),   compareNameTally() ); 
	ofstream fileSummary;
	fileSummary.open(filenameSummary.c_str());

	if (fileSummary.is_open()){

	    fileSummary << "RG\ttotal\ttotal%\tassigned\tassigned%\tunknown\tunknown%\tconflict\tconflict%\twrong\twrong%"<<endl;
	    fileSummary<<dashes<<endl;
	    for(unsigned int i=0;i<toprintVec.size();i++){		
		unsigned int totalForRQ=toprintVec[i].second.assigned+toprintVec[i].second.unknown+toprintVec[i].second.conflict+toprintVec[i].second.wrong;

		fileSummary << toprintVec[i].first << "\t" << totalForRQ << "\t"
			    << 100.0*double(totalForRQ)/double(totalRG) << "%\t" ;

		fileSummary  << toprintVec[i].second.assigned << "\t"
			     << 100.0*double(toprintVec[i].second.assigned)/double(totalForRQ) << "%\t" ;

		fileSummary  << toprintVec[i].second.unknown << "\t"
			     << 100.0*double(toprintVec[i].second.unknown)/double(totalForRQ) << "%\t" ;

		fileSummary <<  toprintVec[i].second.conflict << "\t"
			    << 100.0*double(toprintVec[i].second.conflict)/double(totalForRQ) << "%\t" ;

		fileSummary <<  toprintVec[i].second.wrong << "\t"
			    << 100.0*double(toprintVec[i].second.wrong)/double(totalForRQ) << "%\n" ;
		
		if(toprintVec[i].first != "unknown" )
		    totalAssignRG+=toprintVec[i].second.assigned;
	    }

	    fileSummary<<dashes<<endl;
	    fileSummary<<"ASSIGNED:\t"<< totalAssignRG<<"\t"<<100.0*double(totalAssignRG)/double(totalRG)<<"%"<<endl;
	    fileSummary<<"PROBLEMS:\t"<< (totalRG-totalAssignRG)<<"\t"<<100.0*double(totalRG-totalAssignRG)/double(totalRG)<<"%"<<endl;

	    fileSummary<<"TOTAL:\t"<<totalRG<<"\t100.0%"<<endl;
	}else{
	    cerr << "Unable to print to file "<<filenameSummary<<endl;
	}
	fileSummary.close();
    }




    //Print over-represented sequences in conflict,unknown,wrong
    if(printError){
	vector< pair<string,int> > conflictToPrint( conflictSeq.begin(), conflictSeq.end() ) ;
	vector< pair<string,int> > unknownToPrint(  unknownSeq.begin(),  unknownSeq.end() ) ;
	vector< pair<string,int> > wrongToPrint(    wrongSeq.begin(),    wrongSeq.end() ) ;
     	
	sort (conflictToPrint.begin(),   conflictToPrint.end(),   compareNameRG() ); 
	sort (unknownToPrint.begin(),    unknownToPrint.end(),    compareNameRG() ); 
	sort (wrongToPrint.begin(),      wrongToPrint.end(),      compareNameRG() ); 

	ofstream fileError;
	fileError.open(filenameError.c_str());
	if (fileError.is_open())
	    {
		fileError<<      dashes<<endl<<"Conflict:"<<endl<<dashes<<endl;
		printUnfoundToFile(&conflictToPrint,fileError, printErrorJSONabsent);

		fileError<<endl<<dashes<<endl<<"Unknown:" <<endl<<dashes<<endl;
		printUnfoundToFile(&unknownToPrint,fileError, printErrorJSONabsent);

		fileError<<endl<<dashes<<endl<<"Wrong:"   <<endl<<dashes<<endl;
		printUnfoundToFile(&wrongToPrint,fileError, printErrorJSONabsent);
	    }else{
	    cerr << "Unable to print to file "<<filenameError<<endl;
	}
	fileError.close();
    }

    deallocate();
    if(printError){
	delete trieKnownString;
    }

    if(rgqualFlag)
	rgqualOS.close();

    if(ratioValuesFlag)
	ratioValuesOS.close();
    
    return 0;
}

