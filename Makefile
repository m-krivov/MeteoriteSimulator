# Компиляторы и флаги
CC := gcc
CXX := g++
NVCC := nvcc
CFLAGS := -Wall
CXXFLAGS := -Wall
# -Xptxas=-v для статистики использования ресурсов ядром
NVCCFLAGS := -arch=native --expt-relaxed-constexpr -Xcompiler "-Wall"
LDFLAGS := 
EXEC := app
OBJ_DIR := obj

#INCLUDE_FLAGS := -I. $(patsubst %,-I%,$(wildcard *))
INCLUDE_DIRS := $(shell find . -mindepth 1 -maxdepth 1 -type d)
INCLUDE_FLAGS := $(addprefix -I,$(INCLUDE_DIRS)) -I.

# Поиск исходных файлов (исключаем папку сборки)
C_SOURCES := $(shell find . -type f -name '*.c' ! -path "./$(OBJ_DIR)/*")
CPP_SOURCES := $(shell find . -type f -name '*.cpp' ! -path "./$(OBJ_DIR)/*")
CU_SOURCES := $(shell find . -type f -name '*.cu' ! -path "./$(OBJ_DIR)/*")

# Генерация списка объектных файлов
C_OBJS := $(addprefix $(OBJ_DIR)/, $(C_SOURCES:.c=.o))
CPP_OBJS := $(addprefix $(OBJ_DIR)/, $(CPP_SOURCES:.cpp=.o))
CU_OBJS := $(addprefix $(OBJ_DIR)/, $(CU_SOURCES:.cu=.o))
OBJECTS := $(C_OBJS) $(CPP_OBJS) $(CU_OBJS)

# Автоматическое определение необходимости CUDA
HAS_CUDA := $(if $(CU_SOURCES),yes)

# Правила по умолчанию
all: $(EXEC)

$(EXEC): $(OBJECTS)
ifneq ($(HAS_CUDA),)
	@echo "Linking with CUDA support..."
	$(NVCC) $^ -o $@ $(LDFLAGS)
else
	@echo "Linking without CUDA..."
	$(CXX) $^ -o $@ $(LDFLAGS)
endif

# Компиляция C-файлов
$(OBJ_DIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDE_FLAGS) -c $< -o $@

# Компиляция C++-файлов
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $< -o $@

# Компиляция CUDA-файлов
$(OBJ_DIR)/%.o: %.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) $(INCLUDE_FLAGS) -c $< -o $@

# Очистка
clean:
	rm -rf $(OBJ_DIR) $(EXEC)

.PHONY: all clean