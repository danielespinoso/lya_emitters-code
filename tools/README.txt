
***** HOW MOST COMPLEX FUNCTIONS IN THIS FOLDER WORK *****

 #--#  jplus_reader()  (in data_reader.py)  #--#

This function was written to:
 - read a jplus archive as downloaded from the website, using the Raul's pipeline
 - read a sdss archive as downloaded from CasJobs, using Raul's pipeline
 - crossmatch the two catalogues
 - (eventually) substitute the sdss Broad-Band data to the jplus ones (this last task was needed to
	get a "cleaner", "stable" and "reliable broad band photometry for jplus data, in order to
	operate on them with more confidence)
How it works, step by step, in brief:
- as a first step, it uses Raul's pipeline to read or download the input jplus catalogue (with the
	function 'jplus.datasets.fetch_jplus_objects'. This function requires the input catalogue to
	be in the folder '~/jplus_data/'. If this file is not present, it gets downloaded.
- after that, the sdss catalogue is read. This catalogue must be in the same folder as the routine
	data_reader.py  (i.e.  ~/code/datasets/)
- the scope of having loaded the sdss (photometric) catalogue is to substitute the sdss broad bands
	measurements to the corresponding jplus ones. So the next step is a crossmatch between the two
	catalogues, in order to identify the jplus sources with an sdss counterpart. Once this is
	done, jplus photometric BB magnitudes are replaced by sdss ones.


 #--#  mock_reader()  (in data_reader.py)  #--#

The function was written because:
 - the jplus mock is a huge database containing much more information than my analysis need
 - in order to use on mocks data the same functions I wrote for jplus, I needed to re-format mocks data
 - there are two mocks catalogues: with and without lines. Both contain information that I needed. I
	wanted them in a single input
 - the functions written by David to load mocks are slow
This function was written to:
 - read jplus mock's output
 - select only the needed information
 - write an hdf5 output
How it works, step by step, in brief:
 - it constructs a structure "gal" (since mock data are a collection of "structures" containing the
	information of each object).
 - after that, it defines the fields needed to read mocks data (the can be chosen by de-/commenting)
 - loops over all the 512 light cones containing all data and reads both lines-mock and no-lines-mock.
	NOTE: these two catalogues are in the folder "code_path/datasets/mocks".
 - mocks data come without errors on magnitudes. To introduce them, jplus catalogue is read and
	sampled. After this, new fields in mock catalogue are constructed (namely: "mag_errs")
 - not all mock information is needed, so only narrow-band emitters are selected. The narrow band in
	which the excess has to be measured is an input parameter of the mock_reader routine.
 - the selected data are ultimately formatted in a Raul's pipeline fashion and stored to a hdf5 file.
