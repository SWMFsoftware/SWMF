SHELL =/bin/sh

DEFAULT_TARGET : charge_exchange

include ../../../Makefile.conf

CHARGEEXCHANGEOBJ = ChargeExchange.o ChargeExchange_MaherTinsley.o ChargeExchange_LindsayStebbings.o

charge_exchange : $(CHARGEEXCHANGEOBJ)
	ar -src charge_exchange.a $(CHARGEEXCHANGEOBJ)
	@echo "Target 'charge_exchange' is done." 
