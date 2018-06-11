############################### Default C compiler
default: all


############################### GCC compiler
gcc: CC  = gcc 
gcc: all

############################### Intel compiler
intel: CC  = icc 
intel: all


CFLAGS        = 
INCLUDES      = -I/opt/local/include
LIBS          = -L/opt/local/lib -lgsl -lgslcblas

all:: list gcor4 gmatch pint4

gcor4: 
	 $(CC)  -o $@   $@.c  $(CFLAGS) $(INCLUDES) $(LIBS)

gmatch: 
	 $(CC)  -o $@  $@.c $(CFLAGS) $(INCLUDES) $(LIBS)

pint4: 
	 $(CC)  -o $@   $@.c  $(CFLAGS) $(INCLUDES) $(LIBS)

clean:
	rm -f *.o gcor4 gmatch pint4

help:;
		@echo " "
		@echo "USAGE: "
		@echo "make help                   To get this listing"
		@echo "make                        To compile the program in current environment"
		@echo "make clean                  Remove *.o and executable files"
		@echo "make list                   List the default C compiler "
		@echo "make intel                  To compile the program using the Intel C compiler "
		@echo "make gcc                    To compile the program using the gcc compiler "
		@echo "make [intel,gcc] remake     Fresh make, removing any previous *.o and executable files"
		@echo " "

list:;
		@echo " "
		@echo " Building using the C compiler: $(CC)"
		@echo " "

remake:: clean all

