## Local Setup

1. Clone the repository.
2. Install dependencies: `pip install -r requirements.txt`
3. Set up your local secrets:
   - Go to `.streamlit/`
   - Duplicate `secrets.toml.example` and rename the copy to `secrets.toml`
   - Generate a free access token at [huggingface.co/settings/tokens](https://huggingface.co/settings/tokens)
   - Paste your token into `secrets.toml` under `HF_TOKEN`
4. Run the app: `streamlit run app.py`
