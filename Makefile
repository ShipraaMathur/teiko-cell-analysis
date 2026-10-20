setup:
	python -m pip install --upgrade pip
	python -m pip install -r requirements.txt

pipeline:
	python run_pipeline.py

dashboard:
	python -m streamlit run app.py --server.port 8501 --server.address 0.0.0.0 --server.headless true

clean:
	rm -rf cell_counts.db outputs __pycache__ teiko/__pycache__

.PHONY: setup pipeline dashboard clean
