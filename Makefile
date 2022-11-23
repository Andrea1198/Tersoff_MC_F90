dir=src


all : 
	(cd $(dir) ; make)

clean:
	(cd $(dir) ; make clean)

run:
	(cd $(dir) ; make run)