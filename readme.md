## Local Setup

1. Clone the repository.
2. Install dependencies: `pip install -r requirements.txt`
3. Set up your local secrets:
   - Go to `.streamlit/`
   - Duplicate `secrets.toml.example` and rename the copy to `secrets.toml`
   - Generate a free access token at [huggingface.co/settings/tokens](https://huggingface.co/settings/tokens).
        Under Token type, make sure you select Fine-grained.
        Scroll down to the permissions section, look for Inference, and check the box for Make calls to the serverless Inference API.
   - Paste your token into `secrets.toml` under `HF_TOKEN`
4. Run the app: `streamlit run app.py`
