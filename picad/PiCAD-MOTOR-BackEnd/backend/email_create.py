from pymongo import MongoClient
client = MongoClient('mongodb://localhost:27017/')
db_whitelist = client['whitelist']

# Create the whitelist collection
whitelist_collection = db_whitelist['whitelist_collection']

# List of email IDs to be added to the whitelist
emails_pilabz = [
   "jitin.g@pilabz.in",
  "beeulah.m@pilabz.in",
  "inbaraj.k@pilabz.in",
  "antony.j@pilabz.in",
   "aravindhan.a@pilabz.in",
   "balamuralikrishna.d@pilabz.in",
   "senthilkumar.k@pilabz.in",
   "yuvaraja.t@pilabz.in",
   "ponclinton.a@pilabz.in",
   "oscaranandh.s@pilabz.in",
   "niraipandi.p@pilabz.in",
   "vigneshwar.vw@pilabz.in",
   "supreeth.ss@pilabz.in",
   "venkatesh.sankar@pilabz.in",
   "sharikh.m@pilabz.in",
   "janani.ra@pilabz.in",
   "jagannath.vm@pilabz.in",
   "albin.saji@pilabz.in",
   "manirathnam.ts@pilabz.in",
   "aswanth.mohanan@pilabz.in",
   "pradeep.ak@pilabz.in",
   "balaji.tkk@pilabz.in",
   "mustafa.alam@pilabz.in",
   "chinmaya.vijay@pilabz.in",
   "tangirala.s@pilabz.in",
   "kannan@pilabz.in",
   "hemamalini.j@pilabz.in",
   "guruvignesh@pilabz.in"

]

# Prepare the email documents for insertion
email_docs = [{"email": email} for email in emails_pilabz]

# Insert the email documents into the whitelist collection
whitelist_collection.insert_many(email_docs)