# ------------------------------------------------------------------
# Makefile - CHP TP3 : 3D sampling & confidence intervals
#
# make          compile tout en -O2 (mcPI + bench_O2)
# make bench    compile aussi bench_O0 (sans optimisation)
# make clean    supprime les .o et les executables
# ------------------------------------------------------------------

CC      = gcc
CFLAGS  = -O2 -Wall
CFLAGS0 = -Wall

all: mcPI bench_O2 mc3d

# --- mcPI : exemple 2D fourni (estimation de PI, OpenMP) ----------
mcPI: OMP.c
	$(CC) $(CFLAGS) -fopenmp OMP.c -o mcPI -lm

# --- mc3d : Monte Carlo 3D + intervalles de confiance (OpenMP) ----
mc3d: OMP3D.c mrg32k3a_ctx.c mt19937ar_ctx.c MRG.h mrg32k3a_ctx.h mt19937ar_ctx.h
	$(CC) $(CFLAGS) -fopenmp OMP3D.c mrg32k3a_ctx.c mt19937ar_ctx.c -o mc3d -lm

# --- benchmark des generateurs (MT19937 + MRG de L'Ecuyer) --------
MRG_OBJS = benchmark.o MRG.o mt19937ar.o

bench_O2: $(MRG_OBJS)
	$(CC) $(MRG_OBJS) -o bench_O2

# sans optimisation : recompile depuis les sources, PAS avec les .o -O2
bench_O0:
	$(CC) $(CFLAGS0) benchmark.c MRG.c mt19937ar.c -o bench_O0

benchmark.o: benchmark.c MRG.h mt19937ar.h
	$(CC) $(CFLAGS) -c benchmark.c

MRG.o: MRG.c MRG.h
	$(CC) $(CFLAGS) -c MRG.c

mt19937ar.o: mt19937ar.c mt19937ar.h
	$(CC) $(CFLAGS) -c mt19937ar.c

bench: bench_O0 bench_O2

# --- TODO TP : Monte Carlo 3D + intervalles de confiance ----------
# (implemented: see the mc3d target above)

clean:
	rm -f *.o mcPI bench_O0 bench_O2 mc3d

.PHONY: all bench clean
