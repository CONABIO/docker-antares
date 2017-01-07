# latest version
# execute source /data2/opt/intel/bin/compilervars.sh intel64 before compilation

CC      = gcc
CFLAGS  = -lm -O3  -L/usr/lib64 -lstdc++ -I/usr/include/eigen3/  -I/usr/include/gdal/
LDFLAGS = -lm -fPIC -L/usr/local/lib/ -lgdal -L/usr/lib64 -lstdc++ -ltbb -I/usr/include/eigen3/ 

all: Correlation

correlation: correlation.o
	$(CC) -o $@ $^ $(LDFLAGS)

correlation.o: correlation.cpp 
	$(CC) -c $(CFLAGS) $<

clean:
	rm *.o correlation
