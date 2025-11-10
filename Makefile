PY=python3

all: tables plots

plots:
	$(PY) -m experiments.make_plots

tables:
	$(PY) -m experiments.save_tables

paper:
	pdflatex paper.tex && pdflatex paper.tex

clean:
	rm -f *.aux *.log *.out *.toc *.pdf figures/*.png figures/*.csv figures/*.tex
