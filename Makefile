OBJECTS =  coul.o sle.o  decay.o frag.o random.o array.o loss.o tower.o frame.o tele.o mScat.o fragment.o plf.o pixels.o inverse.o ring.o ttt.o matrix.o
ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
CFLAGS= -c -O2 -std=c++11 -I$(shell root-config --incdir) -g
COMPILER= c++
LINKOPTION = $(shell root-config --libs)


sim: sim.o $(OBJECTS)
	$(COMPILER) -o sim sim.o $(OBJECTS) $(LINKOPTION) 


sim8Be: sim8Be.o $(OBJECTS)
	$(COMPILER) -o sim8Be sim8Be.o $(OBJECTS) $(LINKOPTION) 


2dC: 2dC.o $(OBJECTS)
	$(COMPILER) -o 2dC 2dC.o $(OBJECTS) $(LINKOPTION)

jwkb: jwkb.o $(OBJECTS)
	$(COMPILER) -o jwkb jwkb.o $(OBJECTS) $(LINKOPTION)

b8_6: b8_6.o $(OBJECTS)
	$(COMPILER) -o b8_6 b8_6.o $(OBJECTS) $(LINKOPTION)

be8_22: be8_22.o $(OBJECTS)
	$(COMPILER) -o be8_22 be8_22.o $(OBJECTS) $(LINKOPTION)

b9gs: b9gs.o $(OBJECTS)
	$(COMPILER) -o b9gs b9gs.o $(OBJECTS) $(LINKOPTION)
 
b9_11: b9_11.o $(OBJECTS)
	$(COMPILER) -o b9_11 b9_11.o $(OBJECTS) $(LINKOPTION)
 
$(ALLOBJECTS): %.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@


clean:
	rm -f *.o

