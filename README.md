# Deploy-safe EasyPanel bundle

This bundle includes:
- app.R (fixed v3)
- Dockerfile
- docker-compose.yml
- .dockerignore

## Important
The text you pasted as the “EasyPanel issue” was the app source code, not a build/runtime log.
So this bundle includes safer deployment settings:
- published port 3838:3838
- longer healthcheck start period
- retries increased
- full package list in Dockerfile

## EasyPanel
Deploy from this folder as a Docker Compose project.

## Local Docker
```bash
docker compose down
docker compose build --no-cache
docker compose up -d
```

Then open:
http://YOUR_SERVER_IP:3838

## If it still fails
Check:
```bash
docker logs <container_id>
curl http://127.0.0.1:3838
```
and paste the actual log output.
