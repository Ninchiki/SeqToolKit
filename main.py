from bio_seq import bio_seq

test_DNA = bio_seq()

print(test_DNA.get_seq_info())

DNA = bio_seq("ATCG", 'DNA', 'Test Label')

print(test_DNA.get_seq_info())
test_DNA.gen_random_seq(40, 'DNA')


print(test_DNA.gen_read_frames())

for rf in test_DNA.gen_read_frames():
    print(rf)

print(test_DNA.proteins_from_rf([
    'G', 'M', 'S', 'V', 'G', 'T', 'R', 'R', 'Q', '_', 'A', 'R']))

print(test_DNA.all_proteins_from_RF())










