all: .pclust .pannot

.mpich: mpich.Dockerfile
	docker build --file $< -t darcyabjones/mpich .
	touch $@

.hhsuite: hhsuite.Dockerfile .mpich
	docker build --file $< -t darcyabjones/hhsuite .
	touch $@

.mafft: mafft.Dockerfile
	docker build --file $< -t darcyabjones/mafft .
	touch $@

.mmseqs: mmseqs.Dockerfile .mpich
	docker build --file $< -t darcyabjones/mmseqs .
	touch $@

.python3: python3.Dockerfile
	docker build --file $< -t darcyabjones/python3 .
	touch $@

.pclust: pclust.Dockerfile .hhsuite .mafft .mmseqs .python3
	docker build --file $< -t darcyabjones/pclust .
	touch $@

.emboss: emboss.Dockerfile
	docker build --file $< -t darcyabjones/emboss .
	touch $@

.sperschneider: sperschneider.Dockerfile .emboss
	docker build --file $< -t darcyabjones/sperschneider .
	touch $@

.ffindex: ffindex.Dockerfile .mpich
	docker build --file $< -t darcyabjones/ffindex .
	touch $@

.signalp3: signalp3.Dockerfile signalp-3.0.Linux.tar.Z
	docker build --file $< -t darcyabjones/signalp3 .
	touch $@

.signalp4: signalp4.Dockerfile signalp-4.1f.Linux.tar.gz
	docker build --file $< -t darcyabjones/signalp4 .
	touch $@

.targetp: targetp.Dockerfile targetp-1.1b.Linux.tar.Z chlorop-1.1.Linux.tar.Z .signalp3
	docker build --file $< -t darcyabjones/targetp .
	touch $@

.tmhmm: tmhmm.Dockerfile tmhmm-2.0c.Linux.tar.gz
	docker build --file $< -t darcyabjones/tmhmm .
	touch $@

.phobius: phobius.Dockerfile phobius101_linux.tar.gz
	docker build --file $< -t darcyabjones/phobius .
	touch $@

.pannot: pannot.Dockerfile .python3 .mpich .ffindex .emboss .sperschneider .signalp3 .signalp4 .tmhmm .targetp .phobius
	docker build --file $< -t darcyabjones/pannot .
	touch $@
