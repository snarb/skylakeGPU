################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include ReleseGPU/subdir.mk
-include Release/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(ASM_DEPS)),)
-include $(ASM_DEPS)
endif
ifneq ($(strip $(S_DEPS)),)
-include $(S_DEPS)
endif
ifneq ($(strip $(S_UPPER_DEPS)),)
-include $(S_UPPER_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: IntegrationTests

dependents:
	-cd /home/brans/repos/deepstack/DeepStackCpp/ReleaseSkylakeGPU && $(MAKE) clean all

# Tool invocations
IntegrationTests: $(OBJS) $(USER_OBJS) /home/brans/repos/deepstack/DeepStackCpp/ReleaseSkylakeGPU/libDeepStack.a
	@echo 'Building target: $@'
	@echo 'Invoking: Intel C++ Linker'
	icpc -L/home/brans/repos/deepstack/DeepStackCpp/ReleaseSkylakeGPU/ -L/usr/local/cuda/lib64/ -L"/opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64_lin/" -L/home/brans/repos/deepstack/ziplib/Bin/ -O3 -o "IntegrationTests" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(ASM_DEPS)$(S_DEPS)$(S_UPPER_DEPS)$(CPP_DEPS)$(C_DEPS) IntegrationTests
	-@echo ' '

.PHONY: all clean dependents
/home/brans/repos/deepstack/DeepStackCpp/ReleaseSkylakeGPU/libDeepStack.a:

-include ../makefile.targets