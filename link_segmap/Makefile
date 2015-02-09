# compilation flags
CC = g++
CFLAGS = -ansi $(SPECIAL_FLAGS)
LIBS = -lcfitsio
SRCPATH = .

########
# target section
########
PROGS = $(wildcard $(SRCPATH)/*.cc)
PROGSOBJECTS = $(PROGS:$(SRCPATH)/%.cc=$(SRCPATH)/%)

default: progs

clean: 
	rm -f $(PROGSOBJECTS)

progs: $(PROGSOBJECTS)

$(SRCPATH)/%: $(SRCPATH)/%.cc
	$(CC) $(CFLAGS) $(SPECIAL_LIBS) $< -o $@ $(LIBS)
