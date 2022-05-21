use std::fmt;
use std::error::Error;

#[derive(Debug)]
pub struct MSAError {
    details: String
}

impl MSAError {
    pub fn new(msg: &str) -> MSAError {
        MSAError{details: msg.to_string()}
    }
}

impl fmt::Display for MSAError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,"{}",self.details)
    }
}

impl Error for MSAError {
    fn description(&self) -> &str {
        &self.details
    }
}
