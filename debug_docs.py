#!/usr/bin/env python3

import asyncio
import sys
from pathlib import Path
sys.path.append('src')

from mcp_server.main import handle_read_resource

async def test_docs():
    try:
        print("Testing documentation resources...")

        # Test Docker docs
        docker_docs = await handle_read_resource('documentation://docker')
        print(f'Docker docs length: {len(docker_docs)}')
        print(f'Docker docs content: {docker_docs[:200] if docker_docs else "EMPTY"}')
        print('---')

        # Test Templates
        templates = await handle_read_resource('templates://spatial-workflows')
        print(f'Templates length: {len(templates)}')
        print(f'Templates content: {templates[:200] if templates else "EMPTY"}')
        print('---')

        # Test Nextflow docs
        nextflow_docs = await handle_read_resource('documentation://nextflow')
        print(f'Nextflow docs length: {len(nextflow_docs)}')
        print(f'Nextflow docs content: {nextflow_docs[:200] if nextflow_docs else "EMPTY"}')

    except Exception as e:
        print(f'Error: {e}')
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    asyncio.run(test_docs())
