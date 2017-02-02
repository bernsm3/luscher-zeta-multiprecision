SHELL := /bin/bash

zeta: ; g++ example.cpp -o zeta -Wall -lgmp -lmpfr -O3

plot: ; ./plot.sh

all: zeta

.PHONY: plot zeta
