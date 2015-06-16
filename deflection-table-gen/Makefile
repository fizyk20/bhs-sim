SOURCES       = gen_table.cpp \
		../gr-engine/engine/dpintegrator.cpp \
		../gr-engine/engine/entity.cpp \
		../gr-engine/engine/geometry.cpp \
		../gr-engine/engine/kerr.cpp \
		../gr-engine/engine/kerr_coords.cpp \
		../gr-engine/engine/numeric.cpp \
		../gr-engine/engine/particle.cpp \
		../gr-engine/engine/rk4integrator.cpp \
		../gr-engine/engine/schw.cpp

gen_table: $(SOURCES)
	g++ -I../gr-engine/engine -o gen_table $(SOURCES)