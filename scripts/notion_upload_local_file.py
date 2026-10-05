#!/usr/bin/env python3
"""Send one Notion-prepared local file upload; capability headers stay in stdin.

Use only a fresh create_file_upload response from the authenticated connector.
This helper sends exactly one multipart POST, does not retry or redirect, and
persists only the returned upload receipt, never the capability headers.
"""
import json
import mimetypes
from pathlib import Path
import sys
import termios
import urllib.request
import uuid


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def main():
    if sys.stdin.isatty():
        terminal = termios.tcgetattr(sys.stdin.fileno())
        terminal[3] &= ~(termios.ECHO | termios.ICANON)
        terminal[6][termios.VMIN] = 1
        terminal[6][termios.VTIME] = 0
        termios.tcsetattr(sys.stdin.fileno(), termios.TCSANOW, terminal)
        print('Upload input ready (terminal echo disabled).', flush=True)
    request = json.loads(sys.stdin.readline())
    prepared = request['prepared']
    source = Path(request['source']).resolve()
    receipt = Path(request['receipt']).resolve()
    assert prepared['upload_method'] == 'POST'
    assert prepared['upload_url'].startswith('https://')
    assert source.stat().st_size <= 20 * 1024 * 1024
    filename = prepared['filename']
    assert '\r' not in filename and '\n' not in filename and '"' not in filename
    boundary = 'notion-file-' + uuid.uuid4().hex
    field = prepared.get('upload_form_field', 'file')
    content_type = prepared.get('content_type') or mimetypes.guess_type(filename)[0] or 'application/octet-stream'
    preamble = (f'--{boundary}\r\nContent-Disposition: form-data; name="{field}"; filename="{filename}"\r\nContent-Type: {content_type}\r\n\r\n').encode()
    payload = preamble + source.read_bytes() + f'\r\n--{boundary}--\r\n'.encode()
    headers = dict(prepared['upload_headers'])
    headers['Content-Type'] = f'multipart/form-data; boundary={boundary}'
    req = urllib.request.Request(prepared['upload_url'], data=payload, headers=headers, method='POST')
    opener = urllib.request.build_opener(NoRedirect())
    with opener.open(req, timeout=120) as response:
        result = json.loads(response.read())
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'receipt':str(receipt), 'status':result.get('status'), 'keys':list(result)}))


if __name__ == '__main__':
    main()
