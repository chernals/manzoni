#
# Manzoni
# Global Makefile
#

#
# Define phony target to tell make this is not a file target
.PHONY: manzoni manzonid clean cleanall

#
# Define the manzoni code directory
MANZONI_DIR = manzoni

#
# Make the manzoni complete build
manzoni:
	$(MAKE) -C $(MANZONI_DIR)

manzonid:
	$(MAKE) -C $(MANZONI_DIR)

#
# Install the manzoni general structure
install: manzoni
	mkdir -p ../simulations
	mkdir -p ../simulations/campaigns
	mkdir -p ../simulations/jobs
	chmod +x manzoni/manzoni
	chmod +x pietro/Pietro.rb
	ln -s ../src/manzoni/manzoni ../simulations/manzoni
	ln -s ../src/pietro/Pietro.rb ../simulations/pietro

#
# Set up the examples
examples: manzoni
	ln -s ../src/manzoni/manzoni ../examples/manzoni
	ln -s ../src/pietro/Pietro.rb ../examples/pietro

#
# Other targets
clean:
	$(MAKE) -C $(MANZONI_DIR) clean

cleanall:
	$(MAKE) -C $(MANZONI_DIR) cleanall
