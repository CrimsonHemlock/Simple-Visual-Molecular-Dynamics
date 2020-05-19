#http://clang.llvm.org/docs/AddressSanitizer.html

CXXFLAGS= -O -g -Wextra -Wall -Wshadow -fsanitize=address -Wno-unused-variable -I/opt/local/include -I/opt/X11/include

all: start

start:start.o graphics.o Graphics.h
	$(CXX) $(CXXFLAGS) graphics.o start.o -o $@  -lcairo  -L/opt/X11/lib  -lX11 -lm 


hs.o: hs.cc Graphics.h
start.o: start.cc Graphics.h
graphics.o: Graphics.h graphics.cc

clean:
	$(RM) *.o hs start
clobber: clean
	$(RM)  *.dat projet.tar start hs


tar:
	mkdir MD
	cp Makefile graphics.cc Graphics.h start.cc MD/
	tar cvfz MDprojet.tar MD 
	scp MDprojet.tar acm@turner.pct.espci.fr:public_html/psa135/
	$(RM) -r MD

hs:hs.o graphics.o Graphics.h
	$(CXX) $(CXXFLAGS) graphics.o hs.o -o $@  -L/opt/local/lib -lcairo  -L/opt/X11/lib  -lX11 -lm 
hstar:
	mkdir Attract
	cp Makefile *.cc *.h Attract
	tar cvfz At.tar.gz Attract
