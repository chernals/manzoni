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
	$(MAKE) manzoni -C $(MANZONI_DIR)

manzonid:
	$(MAKE) manzonid -C $(MANZONI_DIR)

#
# Install the manzoni general structure
install: manzoni manzonid
	mkdir -p ../simulations
	mkdir -p ../simulations/campaigns
	mkdir -p ../simulations/jobs
	chmod +x manzoni/manzoni
	chmod +x manzoni/manzonid
	chmod +x pietro/Pietro.rb
	ln -s ../src/manzoni/manzoni ../simulations/manzoni
	ln -s ../src/manzoni/manzonid ../simulations/manzonid
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
