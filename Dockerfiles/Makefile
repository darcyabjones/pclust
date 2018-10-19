docker: .hhblits .sperschneider .signalp3 .signalp4 .targetp .tmhmm .phobius .dssp

.hhblits: hhblits.Dockerfile
	docker build --file hhblits.Dockerfile -t pclust/hhblits .
	touch .hhblits

.sperschneider: sperschneider.Dockerfile
	docker build --file sperschneider.Dockerfile -t pclust/sperschneider .
	touch .sperschneider

.signalp3: signalp3.Dockerfile signalp-3.0.Linux.tar.Z
	docker build --file signalp3.Dockerfile -t pclust/signalp3 .
	touch .signalp3

.signalp4: signalp4.Dockerfile signalp-4.1f.Linux.tar.gz
	docker build --file signalp4.Dockerfile -t pclust/signalp4 .
	touch .signalp4

.targetp: targetp.Dockerfile targetp-1.1b.Linux.tar.Z signalp-3.0.Linux.tar.Z chlorop-1.1.Linux.tar.Z
	docker build --file targetp.Dockerfile -t pclust/targetp .
	touch .targetp

.tmhmm: tmhmm.Dockerfile tmhmm-2.0c.Linux.tar.gz
	docker build --file tmhmm.Dockerfile -t pclust/tmhmm .
	touch .tmhmm

.phobius: phobius.Dockerfile phobius101_linux.tar.gz
	docker build --file phobius.Dockerfile -t pclust/phobius .
	touch .phobius

.dssp: dssp.Dockerfile
	docker build --file dssp.Dockerfile -t pclust/xssp .
	touch .dssp

.blast: blast.Dockerfile
	docker build --file blast.Dockerfile -t pclust/blast .
	touch .blast

.psipred: psipred.Dockerfile .blast .hhblits .dssp
	docker build --file psipred.Dockerfile -t pclust/psipred .
