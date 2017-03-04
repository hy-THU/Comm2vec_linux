CC=g++
CFLAGS=-c -Wall -fpermissive
LDFLAGS= -pthread
SOURCES=C2V_Main.cpp Config.cpp DataSet.cpp Util.cpp Walks.cpp Comm2vec.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=C2V

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o

