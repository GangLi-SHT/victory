include make.inc

ctqmc: 
	@echo '---------------- SRC_MOD ----------------------'
	(cd $(DIR)/src/SRC_MOD; make)
	@echo
	@echo '---------------- SRC_PA --------------------'
	(cd $(DIR)/src/SRC_PA; make)

clean:
	(cd $(DIR)/src/SRC_MOD; make clean)
	(cd $(DIR)/src/SRC_PA; make clean)
	rm -f lib/lib.a

