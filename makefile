TARGET = bin/ipa

CFLAGS = -Wall -pedantic -std=c99 -g
LDFLAGS = -Wall -pedantic -std=c99 -g

SRCS := $(wildcard src/*.c)
OBJS := $(patsubst src/%.c,obj/%.o,$(SRCS))

$(TARGET): $(OBJS)
	@$(CC) $(LDFLAGS) -o $@ $(OBJS)

obj/%.o: src/%.c
	@$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	$(RM) $(OBJS)
