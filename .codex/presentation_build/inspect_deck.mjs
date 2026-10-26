import fs from 'node:fs/promises';
import path from 'node:path';
import { FileBlob, PresentationFile } from '@oai/artifact-tool';
const root = 'C:/Users/yopi1/Documents/09.Research_Dataset/Hanze_SAFEGO';
const source = path.join(root, '00.Meeting/Progress_meeting_Hanze_SAFEGO.pptx');
const out = path.join(root, '.codex/presentation_build/reference');
await fs.mkdir(out, {recursive:true});
const presentation = await PresentationFile.importPptx(await FileBlob.load(source));
const snapshot = await presentation.inspect({kind:'deck,slide,textbox,shape,image,table,chart,notes,layout', maxChars:90000});
await fs.writeFile(path.join(out,'source_inspect.ndjson'), snapshot.ndjson);
const montage = await presentation.export({format:'webp', montage:true, scale:0.8});
await fs.writeFile(path.join(out,'source_montage.webp'), new Uint8Array(await montage.arrayBuffer()));
const slides = presentation.slides.items;
console.log('slide count', slides.length);
console.log(snapshot.ndjson);
for (let i=0;i<slides.length;i++) {
  const preview = await slides[i].export({format:'png', scale:1});
  await fs.writeFile(path.join(out,`slide-${String(i+1).padStart(2,'0')}.png`), new Uint8Array(await preview.arrayBuffer()));
}