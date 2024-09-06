
#[derive (Default, Debug, serde::Deserialize, serde::Serialize)]
pub struct Record {
    pub id: String,
    pub start_byte: u128,
    pub start_bit: usize,
    pub end_byte: u128,
    pub end_bit: usize
}

#[derive (Default, Debug)]
pub struct Records<'a> {
    pub id: &'a str,
    pub records: Vec<Record>,
    pub bytes: Vec<u8>,
    current: usize,
    chunk_size: usize,
    current_rec: usize
}

impl<'a> Records<'a> {
    pub fn new(id: &'a str, records: Vec<Record>, bytes: Vec<u8>, chunk_size: usize) -> Self {
        Self {
            id: id,
            records: records,
            bytes: bytes,
            current: 0,
            chunk_size: chunk_size,
            current_rec: 0
        }
    }
}

impl<'a> Iterator for Records<'a> {
    //type Item = &'a[u8];
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        //TODO: implement
        if self.current == self.bytes.len() - self.chunk_size {
            self.current = 0;
            return None;
        }
        let chunk = self.bytes[self.current..self.current+self.chunk_size].to_vec();
        self.current += 1;
        if self.current as u128 > self.records[self.current_rec].end_byte {
            self.current_rec += 1;
        }
        Some(chunk)
    }

}