# Docker Build and Publish to DockerHub via GitHub Actions

This guide explains how to set up GitHub Actions to automatically build and publish your Docker image to DockerHub.

## Prerequisites

1. A GitHub repository (already set up)
2. A DockerHub account (free at https://hub.docker.com)
3. DockerHub personal access token

## Setup Steps

### 1. Create DockerHub Credentials

1. Log in to [DockerHub](https://hub.docker.com)
2. Navigate to **Account Settings → Security → Personal Access Tokens**
3. Click **Create New Token**
4. Name it (e.g., `github-actions`)
5. Give it **Read & Write** permissions
6. Copy the token (you won't be able to see it again)

### 2. Add GitHub Secrets

In your GitHub repository:

1. Go to **Settings → Secrets and variables → Actions**
2. Click **New repository secret**
3. Add two secrets:
   - **Name:** `DOCKERHUB_USERNAME`
     **Value:** Your DockerHub username
   - **Name:** `DOCKERHUB_TOKEN`
     **Value:** The personal access token you just created

### 3. Workflow Behavior

The workflow file (`.github/workflows/docker-build-publish.yml`) will:

#### On `push` to `main` branch:
- Build the Docker image
- Push to DockerHub with the `latest` tag and commit SHA tag
- Example tags: `davidgllund/plaspart:latest`, `davidgllund/plaspart:abc1234`

#### On `push` of version tags (e.g., `v1.0.0`):
- Build and push with semantic version tags
- Example tags:
  - `davidgllund/plaspart:1.0.0`
  - `davidgllund/plaspart:1.0`
  - `davidgllund/plaspart:latest` (automatic with semver pattern)

#### On Pull Requests:
- Builds the image but does **not** push to DockerHub
- Useful for testing builds without publishing

### 4. Image Naming

The images will be available at:
```
docker.io/your-dockerhub-username/plaspart:tag
```

Examples:
```bash
docker pull docker.io/your-username/plaspart:latest
docker pull docker.io/your-username/plaspart:v1.0.0
```

Or without the registry prefix:
```bash
docker pull your-username/plaspart:latest
```

## Testing the Workflow

1. Push a commit to the `main` branch:
   ```bash
   git push origin main
   ```

2. Or create and push a version tag:
   ```bash
   git tag v1.0.0
   git push origin v1.0.0
   ```

3. Go to your GitHub repository → **Actions** tab
4. Watch the workflow run in real-time
5. Once complete, check your [DockerHub repositories](https://hub.docker.com/repositories)

## Customization

### Change the image name
Edit the `IMAGE_NAME` variable in the workflow file:
```yaml
IMAGE_NAME: my-custom-name
```

### Change the Dockerfile path
Edit the `file` parameter in the build step (currently set to `./app/Dockerfile`)

### Change trigger events
Modify the `on:` section to trigger on different events:
- `workflow_dispatch` - manual trigger via GitHub UI
- `schedule` - cron-based builds
- `release` - on GitHub releases

## Troubleshooting

### Workflow fails with authentication error
- Verify `DOCKERHUB_USERNAME` and `DOCKERHUB_TOKEN` are set correctly in GitHub Secrets
- Ensure the token is active and hasn't expired

### Build fails
- Check the workflow logs in GitHub Actions
- Test building locally: `docker build -f app/Dockerfile -t plaspart:test .`
- Ensure all required files are in the correct paths

### Image not appearing on DockerHub
- Check if the workflow ran successfully (green checkmark)
- Verify you're looking at the correct DockerHub account
- Check that the repository is set to public if needed

## Advanced Options

### Push to multiple registries
You can extend the workflow to push to both DockerHub and GitHub Container Registry (GHCR).

### Build for multiple architectures
Add this to the `build-push-action` step:
```yaml
platforms: linux/amd64,linux/arm64
```

### Skip DockerHub on specific commits
Add `[skip docker]` to your commit message if using a conditional step.

## Next Steps

- Configure repository-level build rules in DockerHub to pull from GitHub
- Set up automated scanning for vulnerabilities
- Monitor image size and optimize the Dockerfile if needed
