import os, sys

def readNextFourLines(fp):

  lines = []
  for n in range(4):
    line = fp.readline()
    if not line:
      return []
    lines.append(line)

  return lines

def writeLines(fp, lines, barcodeLen):

  n = barcodeLen + 1
  #n = barcodeLen
  fp.write(lines[0])
  fp.write(lines[1][n:])
  fp.write(lines[2])
  fp.write(lines[3][n:])

def main(inputFile1, inputFile2, outputDirectory, barcodeLen=3):

  if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

  fpin1 = open(inputFile1, 'rU')
  fpin2 = open(inputFile2, 'rU')

  barcodeCount = {}
  barcodeFp1 = {}
  barcodeFp2 = {}
  mismatched = 0
  matched = 0
  Ncount1 = 0
  Ncount2 = 0

  lines1 = readNextFourLines(fpin1)
  while lines1:
    if lines1[0][0] != '@':
      raise Exception('lines1[0] does not start with @: "%s"' % lines1[0].strip)
    lines2 = readNextFourLines(fpin2)
    if lines2[0][0] != '@':
      raise Exception('lines2[0] does not start with @: "%s"' % lines2[0].strip)

    if lines1[1][barcodeLen] == 'N':
      Ncount1 += 1

    if lines2[1][barcodeLen] == 'N':
      Ncount2 += 1

    barcode1 = lines1[1][:barcodeLen]
    barcode2 = lines2[1][:barcodeLen]
    if barcode1 != barcode2:
      mismatched += 1
    else:
      matched += 1
      if barcode1 not in barcodeFp1:
        fileName1 = os.path.join(outputDirectory, '%s_%s.fq' % (os.path.basename(inputFile1)[:-3], barcode1))
        fileName2 = os.path.join(outputDirectory, '%s_%s.fq' % (os.path.basename(inputFile2)[:-3], barcode2))
        fpout1 = open(fileName1, 'w')
        fpout2 = open(fileName2, 'w')
        barcodeFp1[barcode1] = fpout1
        barcodeFp2[barcode2] = fpout2
        barcodeCount[barcode1] = 0

      fpout1 = barcodeFp1[barcode1]
      fpout2 = barcodeFp2[barcode2]
      barcodeCount[barcode1] += 1

      writeLines(fpout1, lines1, barcodeLen)
      writeLines(fpout2, lines2, barcodeLen)
        
    lines1 = readNextFourLines(fpin1)

  fpin1.close()
  fpin2.close()

  for fpout1 in barcodeFp1.values():
    fpout1.close()

  for fpout2 in barcodeFp2.values():
    fpout2.close()

  total = matched + mismatched
  print('total = %d' % total)
  total = float(total)
  print('matched = %d (%.1f%%)' % (matched, 100*matched/total))
  print('mismatched = %d (%.1f%%)' % (mismatched, 100*mismatched/total))
  print('N at fourth character in first file = %d (%.1f%%)' % (Ncount1, 100*Ncount1/total))
  print('N at fourth character in second file = %d (%.1f%%)' % (Ncount2, 100*Ncount2/total))

  for barcode in sorted(barcodeCount):
    count = barcodeCount[barcode]
    print('barcode "%s" = %d (%.1f%%)' % (barcode, count, 100*count/total))

if __name__ == '__main__':
  
  in_files = sys.argv[1:]
  in_file_1, in_file2 = in_files[:2]
  out_dir = os.path.dirname(os.path.abspath(in_file_1))
  
  main(in_file_1, in_file2, out_dir)
