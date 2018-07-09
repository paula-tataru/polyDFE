############################################################
# for non-default GSL installation
############################################################
USR_INC := -IC:/cygwin64/usr/local/include/
USR_LIB := -LC:/cygwin64/usr/local/lib/
############################################################

RM := rm -rf

LIBS := -lm -lgsl -lgslcblas

OPT := -O2 -g3 -Wall -c -msse2 -mfpmath=sse -frounding-math -fmessage-length=0 -MMD -MP

C_SRCS += \
../src/basinhopping.c \
../src/likelihood.c \
../src/localoptim.c \
../src/parse.c \
../src/run.c \
../src/simulate.c \
../src/transform.c

OBJS += \
./src/basinhopping.o \
./src/likelihood.o \
./src/localoptim.o \
./src/parse.o \
./src/run.o \
./src/simulate.o \
./src/transform.o

C_DEPS += \
./src/basinhopping.d \
./src/likelihood.d \
./src/localoptim.d \
./src/parse.d \
./src/run.d \
./src/simulate.d \
./src/transform.d


%.o: %.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc $(USR_INC) $(OPT) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


all: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc -o "polyDFE" $(USR_LIB) $(OBJS) $(LIBS)
	@echo 'Finished building target: $@'


clean:
	-$(RM) $(EXECUTABLES) $(OBJS) $(C_DEPS) polyDFE
	-@echo ' '
