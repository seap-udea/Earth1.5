
#LIBDIR=/Users/jzuluaga/usr
LIBDIR=/usr/local
LIBGSL=$(LIBDIR)
LIBCONFIG=util

CPP=g++
# LFLAGS=-lm -L$(LIBDIR)/lib -lgsl -lgslcblas -L$(LIBCONFIG)/lib/ -lconfig++
LFLAGS=-lm $(LIBDIR)/lib/libgsl.a $(LIBDIR)/lib/libgslcblas.a $(LIBCONFIG)/lib/libconfig++.a
CFLAGS=-c -I. -I$(LIBGSL)/include -I$(LIBCONFIG)/include

%.out:%.o
	$(CPP) $^ $(LFLAGS) -o $@

%.o:%.cpp
	$(CPP) $(CFLAGS) $^ -o $@

clean:
	rm *.out

commit:
	@echo "Commiting changes..."
	@git commit -am "Commit"
	@git push origin master

pull:
	@echo "Pulling from repository..."
	@git reset --hard HEAD	
	@git pull
