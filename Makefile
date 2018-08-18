
#replace with the absolute path to the folder /nauty26r7
#download from http://pallini.di.uniroma1.it/#howtogetit and install

NAUTYHOME := $(HOME)/.opt/nauty

CXX := clang++
CXXFLAGS := -std=c++14

HEADER_DIR := h

NAUTY_LIB := $(NAUTYHOME)/nauty.a

prog:
	$(CXX) $(CXXFLAGS) -O3 -g 2L_enum.cpp -o 2L_enum $(NAUTY_LIB) -I$(NAUTYHOME)