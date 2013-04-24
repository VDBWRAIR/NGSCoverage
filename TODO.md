TODO LIST
=========

* Since I am quite new at matplotlib I have been quite unsuccessful at decoupling the pyplt calls. That is, it would be
  much nicer if each gap/low coverage region was an object of its own that would graph its own line. Aka, each region
  would be a pyplot.figure(?) object or maybe an Axes object, I'm not 100% sure how everything interacts in Matplotlib.
* Would be useful in the future to somehow be able to distinquish segments from an entire genome and graph each separtely
  on one graph each in its own subplot. Right now, you can only achieve something similar by splitting out the segments
  manually(or through scripting) and then graphing each separately and then joining the images manually
* Model everything after a more generic alignment. That is, instead of focusing on GsMapper projects, it would be a lot
  better if the BAM pileups were used. I feel that may be a duplication of effort as I think there may be a project
  already out there that does this.
