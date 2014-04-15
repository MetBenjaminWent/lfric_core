##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

DOXYGEN ?= doxygen

.PHONY: api
api:
	$(DOXYGEN) Doxyfile

REPOSITORY_URL = https://puma.nerc.ac.uk/svn/GungHo_svn/Dynamo

.PHONY: gource
gource: gource/head.txt dynamo-gource.mp4

dynamo-gource.mp4: gource/log.xml gource/captions.txt
	$(info *** The 'gource' visualisation tool takes frames straight from rendered output.)
	$(info *** If you obscure the window during rendering frames will be incomplete.)
	gource -1280x720 -r 30 --multi-sampling --title Dynamo \
	       --user-image-dir avatars --caption-file gource/captions.txt \
	       --hide mouse,progress,usernames,filenames,dirnames -o - \
               gource/log.xml \
	       | ffmpeg -y -r 30 -f image2pipe -vcodec ppm -i - -c:v mpeg4 \
	       -s 1280x720 -r 30 -pix_fmt yuv420p -g 15 -b:v 2000k \
	       -maxrate 4000k -minrate 1000k -bufsize 5000k -packetsize 4k \
	       -muxrate 5000k $@

gource/log.xml: gource/head.txt
	svn log -r 1:HEAD --xml --verbose --quiet $(REPOSITORY_URL) > $@

gource/captions.txt: gource/head.txt
	tools/encaptor $(REPOSITORY_URL)/trunk $@

gource/head.txt: FORCE
	tools/headaliser $(REPOSITORY_URL) gource/head.txt

.PHONY: codeswarm
codeswarm: code_swarm/head.txt dynamo-codeswarm.mp4

dynamo-codeswarm.mp4: code_swarm/log.xml
	code_swarm -c code_swarm/dynamo.config

code_swarm/log.xml: code_swarm/log.txt
	convert_logs.py --svn-log=$< --output=$@

code_swarm/log.txt:
	svn log -v $(REPOSITORY_URL) > $@

code_swarm/head.txt: FORCE
	tools/headaliser $(REPOSITORY_URL) code_swarm/head.txt

clean: FORCE
	-rm -rf docs/api/*
	-rm gource/head.txt gource/captions.txt gource/log.xml dynamo-gource.mp4

FORCE:
