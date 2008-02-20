
CXX = g++

DEBUG =
#-g

###############################################################################

IFLAGS = $(DEBUG) $(INCLUDE)
FLAGS = $(DEBUG) $(LIB)

OBJECTS = eke.o grid.o

###############################################################################

all: eke

clean:
	@rm -rfv eke *.o *.so

distclean: clean
	@rm -rfv *~

%.o: %.cpp *.hpp
	$(CXX) $(IFLAGS) -c $<

eke: $(OBJECTS) 
	$(CXX) -o eke $(OBJECTS)

