all: mmseqs2 sperschneider signalp3 signalp4 targetp tmhmm phobius

hhblits: hhblits.Dockerfile
	docker build --file hhblits.Dockerfile -t pclust/hhblits .

sperschneider: sperschneider.Dockerfile
	docker build --file sperschneider.Dockerfile -t pclust/sperschneider .

signalp3: signalp3.Dockerfile signalp-3.0.Linux.tar.Z
	docker build --file signalp3.Dockerfile -t pclust/signalp3 .

signalp4: signalp4.Dockerfile signalp-4.1f.Linux.tar.gz
	docker build --file signalp4.Dockerfile -t pclust/signalp4 .

targetp: targetp.Dockerfile targetp-1.1b.Linux.tar.Z signalp-3.0.Linux.tar.Z chlorop-1.1.Linux.tar.Z
	docker build --file targetp.Dockerfile -t pclust/targetp .

tmhmm: tmhmm.Dockerfile tmhmm-2.0c.Linux.tar.gz
	docker build --file tmhmm.Dockerfile -t pclust/tmhmm .

phobius: phobius.Dockerfile phobius101_linux.tar.gz
	docker build --file phobius.Dockerfile -t pclust/phobius .

