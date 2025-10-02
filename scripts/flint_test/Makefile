CC = gcc
# Add common optimizations and OpenMP
CFLAGS = -Wall -O3 -fopenmp
# Ensure headers and libs from /usr/local are considered (common FLINT install prefix)
CFLAGS += -I/usr/local/include
LDFLAGS = -L/usr/local/lib -lcrypto -lm -lflint -lgmp
# Embed rpath so the runtime linker can find libflint without LD_LIBRARY_PATH
LDFLAGS += -Wl,-rpath,/usr/local/lib

TARGET = paillier_bench
SRCS = paillier_bench.c
OBJS = $(SRCS:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
