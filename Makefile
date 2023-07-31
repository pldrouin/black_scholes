CNAMES := black_scholes.c root_finder.c

LCOBJ := $(CNAMES:.c=.o)
LCODEP := $(LCOBJ:.o=.d)

LIBNAME:= black_scholes
LIB	:= lib$(LIBNAME).so

all:	$(LIB)

$(LIB): $(LCOBJ)
	$(CC) $(CFLAGS) -fPIC -shared -o $@ $^ $(LDFLAGS) -lm

$(LCODEP): %.d: %.c %.h
	@echo "Generating dependency file $@"
	@set -e; rm -f $@
	@$(CC) -M $(CFLAGS) $< > $@.tmp
	@sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.tmp > $@
	@rm -f $@.tmp

include $(LCODEP)

$(LCOBJ): %.o: %.c %.h
	$(CC) $(CFLAGS) -fPIC -c -o $@ $<

clean:
	rm -rf $(LCOBJ) $(LCODEP)

clear: clean
	rm -rf $(LIB)
