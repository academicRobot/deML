
BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)
SIMPLEJSON= $(realpath ../SimpleJSON/)



CXXFLAGS+=-std=c++11 -O3 -Wall -I${BAMTOOLS}/src/ -I${LIBGAB} -I${LIBGAB}/gzstream/  -I${SIMPLEJSON}/src/
LDLIBS+=/usr/lib/x86_64-linux-gnu/libbamtools.a ${LIBGAB}/libgab.a ${LIBGAB}/gzstream/libgzstream.a  -lm -lz

CXX=g++ #-g -pg

all:    deML

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

deML: deML.o PrefixTree.o RGAssign.o ${LIBGAB}/PutProgramInHeader.o ${SIMPLEJSON}/obj/JSON.o ${SIMPLEJSON}/obj/JSONValue.o ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o 
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

static: deML.o PrefixTree.o RGAssign.o ${LIBGAB}/PutProgramInHeader.o ${SIMPLEJSON}/obj/JSON.o ${SIMPLEJSON}/obj/JSONValue.o ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o 
	${CXX} -static $(LDFLAGS) -o deML $^ $(LDLIBS)



clean:
	rm -f deML *.o

