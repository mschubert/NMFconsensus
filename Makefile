all:
	$(MAKE) -C libnmf/ libnmf.so
	ln -fs libnmf/libnmf.so libnmf.so

clean:
	$(MAKE) -C libnmf/ clean
	rm -f libnmf.so

