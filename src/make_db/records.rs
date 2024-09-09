use std::sync::Arc;
#[derive (Default, Debug, serde::Deserialize, serde::Serialize)]
pub struct Record {
    pub id: String,
    pub start_byte: u128,
    pub start_bit: usize,
    pub end_byte: u128,
    pub end_bit: usize
}

#[derive (Default)]
pub struct VecRecord {
    current: usize,
    pub current_record: Arc<Record>,
    records: Vec<Arc<Record>>
}

impl Iterator for VecRecord {
    type Item = Arc<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current == self.records.len() {
            return None;
        }
        self.current_record = Arc::clone(&self.records[self.current]);
        self.current += 1;
        Some(Arc::clone(&self.current_record))
    }
}

impl VecRecord {
    pub fn new(records: Vec<Arc<Record>>) -> Self {
        Self {
            current_record : Arc::clone(&records[0]),
            records: records,
            current: 0,
        }
    }

    pub fn push(&mut self, record: Record) {
        self.records.push(Arc::new(record));
    }
}

#[derive (Default, Debug)]
pub struct Records<'a> {
    pub id: &'a str,
    pub records: Vec<Record>,
    pub bytes: Vec<u8>,
    current: usize,
    chunk_size: usize,
    pub current_rec: usize
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

pub struct LazyRecords<'a> {
    pub id: &'a str,
    pub records: Vec<Record>,

}

pub struct SimpleRecord {
    pub seq: Vec<u8>,
    pub id: String,
    pub words: Vec<Vec<u8>>,
    pub k: usize
}