SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g
LIBLOR = $(HOME_LORENE)/Lib/liblorene_g.a $(HOME_LORENE)/Lib/liblorenef77_g.a

.SUFFIXES : .o .C

SRC1 = file_reader.C
OBJ1 = $(SRC1:.C=.o)
EXE1 = $(SRC1:.C=)




$(EXE1): $(OBJ1) $(LIBLOR)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
	cp $(HOME_LORENE)/Visu/OpenDX/visu_scalar_box.net .
	cp $(HOME_LORENE)/Visu/OpenDX/visu_section.net .
	cp $(HOME_LORENE)/Visu/OpenDX/visu_vector3d.net .


# $(LIB_PGPLOT)

all: $(EXE1)


.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

uninstall:
	rm -f $(OBJ1)
	rm -f $(EXE1)
	rm -f visu*.net
